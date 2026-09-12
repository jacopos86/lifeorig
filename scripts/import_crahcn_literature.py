#!/usr/bin/env python3
"""Convert the published CRAHCN LaTeX tables to LIFEORIG reaction files."""

import argparse
import re
from pathlib import Path


NUMBER = re.compile(r"\*?(\d+[a-z]?)\.\s*&")


def table(text, label):
    marker = f"\\label{{{label}}}"
    start = text.index(marker)
    start = text.index("\\\\", start) + 2
    end = text.index("\\end{longtable*}", start)
    return text[start:end]


def records(table_text, numbered=True):
    result = []
    current = None
    for line in table_text.splitlines():
        match = NUMBER.match(line.strip()) if numbered else None
        if match:
            if current:
                result.append(current)
            current = {"id": match.group(1), "lines": [line]}
        elif numbered and current:
            current["lines"].append(line)
        elif not numbered and "\\ce{" in line and line.lstrip().startswith("\\ce{"):
            result.append({"id": str(len(result) + 1), "lines": [line]})
    if current:
        result.append(current)
    return result


def fields(line):
    return [part.strip() for part in line.rstrip(" \\").split("&")]


def ce_from_reaction_column(lines):
    expressions = []
    for line in lines:
        parts = fields(line)
        if len(parts) < 2:
            continue
        match = re.search(r"\\ce\{(.+)\}", parts[1])
        if match:
            expressions.append(match.group(1))
    return " ".join(expressions)


def normalize_species(text):
    replacements = {
        "^1CH2": "CH21",
        "^3CH2": "CH23",
        "^1O": "O1D",
        "^3O": "O",
        "^4N": "N",
        "^2N": "N2D",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"\$_\{\(\\nu\)\}\$", "", text)
    text = re.sub(r"_\{\(\\nu\)\}", "", text)
    text = text.replace("$\\cdots$", "").replace("*", "")
    text = re.sub(r"\^\{[^}]+\}", "", text)
    text = text.replace("{", "").replace("}", "")
    return " ".join(text.split())


def net_equation(expression):
    expression = normalize_species(expression)
    if "<-" in expression and "->" not in expression:
        pieces = [piece.strip() for piece in expression.split("<-")]
        lhs, rhs = pieces[-1], pieces[0]
    else:
        pieces = [piece.strip() for piece in expression.split("->")]
        lhs, rhs = pieces[0], pieces[-1]
    if not lhs or not rhs:
        raise ValueError(f"Cannot reduce reaction pathway: {expression}")
    return f"{lhs} -> {rhs}"


def latex_number(value):
    value = value.strip()
    match = re.search(r"([+-]?\d+(?:\.\d+)?)\$\\times\$10\$\^\{([+-]?\d+)\}\$", value)
    if match:
        return f"{match.group(1)}e{match.group(2)}"
    match = re.search(r"[+-]?\d+(?:\.\d+)?", value)
    return match.group(0) if match else ""


def make_rows(source, table_label, module, prefix, paper, kind):
    parsed = records(table(source, table_label))
    rows = []
    for record in parsed:
        first = fields(record["lines"][0])
        expression = ce_from_reaction_column(record["lines"])
        if not expression:
            continue
        # CRAHCN reaction 74 is split across two source-table rows without
        # an arrow before the final products. Preserve its published net path.
        if prefix == "CR" and table_label == "CRAHCN2" and record["id"] == "74":
            expression = "^2N + CH3 -> H2CN + H"
        if kind == "lindemann":
            k_inf = latex_number(first[2])
            k_zero = latex_number(first[3])
            rate = f"lindemann(k_inf={k_inf},k0={k_zero},T_ref=298K)"
            control = "three_body"
        else:
            alpha = latex_number(first[3])
            beta = latex_number(first[4])
            gamma = latex_number(first[5])
            if not alpha:
                continuation = [fields(line) for line in record["lines"][1:]]
                parameter_rows = [row for row in continuation if len(row) > 5 and latex_number(row[3])]
                values = [
                    f"alpha={latex_number(row[3])},beta={latex_number(row[4])},gamma={latex_number(row[5])}"
                    for row in parameter_rows
                ]
                rate = "piecewise_modified_arrhenius(" + ";".join(values) + ")"
            else:
                rate = f"modified_arrhenius(alpha={alpha},beta={beta},gamma={gamma},T_ref=300K)"
            control = "thermal"
        reaction_id = f"{prefix}{record['id'].upper()}"
        reference = f"{paper}; Table {table_label}, reaction {record['id']}"
        rows.append(
            " | ".join(
                [reaction_id, module, net_equation(expression), control, rate,
                 "atmospheric gas-phase kinetics", reference, "published_network"]
            )
        )
    return rows


def make_hadean_rows(source):
    parsed = records(table(source, "addedtwobodychemistry"), numbered=False)
    rows = []
    for record in parsed:
        first = fields(record["lines"][0])
        expression = re.search(r"\\ce\{(.+)\}", first[0]).group(1)
        alpha = latex_number(first[1])
        beta = latex_number(first[2])
        gamma = latex_number(first[3])
        rate = f"modified_arrhenius(alpha={alpha},beta={beta},gamma={gamma},T_ref=300K)"
        rows.append(" | ".join([
            f"HAD{int(record['id']):03d}", "pearce2022_hadean_extension",
            net_equation(expression), "thermal", rate,
            "Hadean gas-phase extension",
            f"Pearce, He & Horst 2022, DOI 10.1021/acsearthspacechem.2c00138; reaction {record['id']}",
            "published_network",
        ]))
    return rows


def write_network(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + "\n\n" + "\n".join(rows) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crahcn", type=Path, required=True)
    parser.add_argument("--crahcn-o", type=Path, required=True)
    parser.add_argument("--hadean", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    crahcn = args.crahcn.read_text(encoding="utf-8")
    oxygen = args.crahcn_o.read_text(encoding="utf-8")
    hadean = args.hadean.read_text(encoding="utf-8")
    base_paper = "Pearce et al. 2020, DOI 10.3847/1538-4357/abae5b"
    oxygen_paper = "Pearce et al. 2020, DOI 10.1021/acs.jpca.0c06804"
    rows = (
        make_rows(crahcn, "CRAHCN1", "crahcn_lindemann", "CR", base_paper, "lindemann")
        + make_rows(crahcn, "CRAHCN2", "crahcn_arrhenius", "CR", base_paper, "arrhenius")
        + make_rows(oxygen, "CRAHCN1", "crahcn_o_lindemann", "CRO", oxygen_paper, "lindemann")
        + make_rows(oxygen, "CRAHCN2", "crahcn_o_arrhenius", "CRO", oxygen_paper, "arrhenius")
    )
    header = """# Complete CRAHCN plus CRAHCN-O oxygen extension, transcribed from the published supporting tables.
# Columns: ID | MODULE | REACTION | CONTROL | RATE_TEMPLATE | ROLE | REFS | CONFIDENCE
# Multi-step mechanisms are reduced to their net reactants and final products; source files remain unmodified."""
    write_network(args.output / "crahcn_o_pearce2020.txt", header, rows)

    extension = make_hadean_rows(hadean)
    header = """# 53 two-body Hadean-atmosphere reactions from Pearce, He & Horst (2022).
# Columns: ID | MODULE | REACTION | CONTROL | RATE_TEMPLATE | ROLE | REFS | CONFIDENCE
# Multi-step mechanisms are reduced to their net reactants and final products; source files remain unmodified."""
    write_network(args.output / "hadean_nh3_no_extension_pearce2022.txt", header, extension)
    print(f"CRAHCN-O: {len(rows)} reactions")
    print(f"Hadean extension: {len(extension)} reactions")


if __name__ == "__main__":
    main()
