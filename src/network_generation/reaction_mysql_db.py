import os
import pymysql
from collections import Counter
from src.utilities.logging_module import log

def open_reaction_database():
    return pymysql.connect(
        host=os.environ["LIFEORIG_REACTION_DB_HOST"],
        user=os.environ["LIFEORIG_REACTION_DB_USER"],
        password=os.environ["LIFEORIG_REACTION_DB_PASSWORD"],
        database=os.environ["LIFEORIG_REACTION_DB_NAME"],
        autocommit=False,
    )

def reaction_file_is_current(db, reaction_file, file_hash):
    with db.cursor() as cursor:
        cursor.execute(
            """
            SELECT file_hash
            FROM reaction_source_files
            WHERE file_path = %s
            """,
            (str(reaction_file),),
        )
        row = cursor.fetchone()
    return row is not None and row[0] == file_hash

def delete_old_reaction_file_data(db, reaction_file):
    with db.cursor() as cursor:
        cursor.execute(
            """
            DELETE FROM reaction_source_files
            WHERE file_path = %s
            """,
            (str(reaction_file),),
        )

#
#   insert reactions file utilities
#

def insert_reaction_file_data(db, reaction_file, file_hash, parsed_network):
    with db.cursor() as cursor:
        source_file_id = _insert_reaction_source_file(
            cursor, reaction_file, file_hash
        )
        # species id
        species_ids = {}
        for species_name in parsed_network.species:
            species_ids[species_name] = _get_or_insert_species(
                cursor, species_name
            )
        # reactions
        for reaction in parsed_network.reactions:
            reaction_id = _insert_reaction(
                cursor, source_file_id, reaction
            )
            _insert_reaction_species(
                cursor,
                reaction_id,
                reaction.reactants,
                species_ids,
                "reactant",
            )
            _insert_reaction_species(
                cursor,
                reaction_id,
                reaction.products,
                species_ids,
                "product",
            )

def _insert_reaction_source_file(cursor, reaction_file, file_hash):
    cursor.execute(
        """
        INSERT INTO reaction_source_files (environment_model, file_path, file_hash)
        VALUES (%s, %s, %s)
        """,
        ("unknown", str(reaction_file), file_hash),
    )
    return cursor.lastrowid

def _get_or_insert_species(cursor, species_name):
    cursor.execute(
        """
        INSERT INTO species (name)
        VALUES (%s)
        ON DUPLICATE KEY UPDATE id = LAST_INSERT_ID(id)
        """,
        (species_name,),
    )
    return cursor.lastrowid

def _insert_reaction(cursor, source_file_id, reaction):
    cursor.execute(
        """
        INSERT INTO reactions (
            source_file_id,
            source_reaction_id,
            module,
            equation,
            reversible,
            catalyst_or_control,
            rate_template,
            role,
            refs,
            confidence
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """,
        (
            source_file_id,
            reaction.reaction_id,
            reaction.module,
            reaction.equation,
            reaction.reversible,
            reaction.catalyst_or_control,
            reaction.rate_template,
            reaction.role,
            reaction.refs,
            reaction.confidence,
        ),
    )
    return cursor.lastrowid

def _insert_reaction_species(cursor, reaction_id, species_names, species_ids, side):
    for species_name, stoichiometry in Counter(species_names).items():
        cursor.execute(
            """
            INSERT INTO reaction_participants (
                reaction_id,
                species_id,
                side,
                stoichiometry
            )
            VALUES (%s, %s, %s, %s)
            """,
            (
                reaction_id,
                species_ids[species_name],
                side,
                float(stoichiometry),
            ),
        )

def _source_file_filter(source_files):
    if not source_files:
        return "", []
    paths = [str(source_file) for source_file in source_files]
    placeholders = ", ".join(["%s"] * len(paths))
    return f"WHERE rsf.file_path IN ({placeholders})", paths


def get_species_from_reaction_database(db, source_files=None):
    source_filter, source_params = _source_file_filter(source_files)
    with db.cursor() as cursor:
        cursor.execute(
            f"""
            SELECT DISTINCT s.name
            FROM species s
            JOIN reaction_participants rp ON rp.species_id = s.id
            JOIN reactions r ON r.id = rp.reaction_id
            JOIN reaction_source_files rsf ON rsf.id = r.source_file_id
            {source_filter}
            ORDER BY s.name
            """,
            source_params,
        )
        return [row[0] for row in cursor.fetchall()]

#
#   SUMMARY INFO
#

def log_species_summary(db, source_files=None):
    source_filter, source_params = _source_file_filter(source_files)
    with db.cursor() as cursor:
        cursor.execute(
            f"""
            SELECT COUNT(*)
            FROM reaction_source_files rsf
            {source_filter}
            """,
            source_params,
        )
        n_source_files = cursor.fetchone()[0]
        cursor.execute(
            f"""
            SELECT COUNT(*)
            FROM reactions r
            JOIN reaction_source_files rsf ON rsf.id = r.source_file_id
            {source_filter}
            """,
            source_params,
        )
        n_reactions = cursor.fetchone()[0]
        cursor.execute(
            f"""
            SELECT DISTINCT s.name
            FROM species s
            JOIN reaction_participants rp ON rp.species_id = s.id
            JOIN reactions r ON r.id = rp.reaction_id
            JOIN reaction_source_files rsf ON rsf.id = r.source_file_id
            {source_filter}
            ORDER BY s.name
            """,
            source_params,
        )
        species = [row[0] for row in cursor.fetchall()]
    log.info("\t REACTION DATABASE SUMMARY")
    log.info("\t source files : " + str(n_source_files))
    log.info("\t reactions    : " + str(n_reactions))
    log.info("\t species      : " + str(len(species)))
    log.info("\t species list : " + ", ".join(species))
