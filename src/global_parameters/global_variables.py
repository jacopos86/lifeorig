from pathlib import Path
import src

PACKAGE_DIR = Path(src.__path__[0]).resolve()
PROJECT_DIR = PACKAGE_DIR.parent