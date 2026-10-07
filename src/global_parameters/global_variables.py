from pathlib import Path
import src

PACKAGE_DIR = Path(src.__path__[0]).resolve()
PROJECT_DIR = PACKAGE_DIR.parent

# SOURCE REACTION FILES

VPL_ARCHEAN_REACTIONS = (
    Path(__file__).resolve().parents[2]
    / "reference_reactions"
    / "ATMOSPHERIC"
    / "VPL_ATMOS"
    / "PHOTOCHEM"
    / "INPUTFILES"
    / "TEMPLATES"
    / "Archean+haze"
    / "reactions.rx"
).resolve()