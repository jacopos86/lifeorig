import re
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


SRC = Path("doc/ATMOSPHERE/atmosphere_hydrothermal_architecture.tex")
OUT = Path("doc/ATMOSPHERE/atmosphere_hydrothermal_architecture.pdf")


def latex_to_lines(raw: str) -> list[str]:
    body = raw.split(r"\begin{document}", 1)[1].split(r"\end{document}", 1)[0]
    lines = []
    in_code = False
    for line in body.splitlines():
        s = line.rstrip()
        if s.startswith(r"\maketitle"):
            continue
        if s.startswith(r"\begin{lstlisting}"):
            in_code = True
            lines.append("")
            continue
        if s.startswith(r"\end{lstlisting}"):
            in_code = False
            lines.append("")
            continue
        if in_code:
            lines.append("    " + s)
            continue

        s = re.sub(r"\\section\{([^}]*)\}", r"\n\1\n" + "=" * 72, s)
        s = re.sub(r"\\subsection\{([^}]*)\}", r"\n\1\n" + "-" * 48, s)
        for token in (
            r"\begin{itemize}",
            r"\end{itemize}",
            r"\begin{enumerate}",
            r"\end{enumerate}",
            r"\begin{description}[leftmargin=3.7cm]",
            r"\end{description}",
            r"\begin{quote}",
            r"\end{quote}",
        ):
            s = s.replace(token, "")
        s = re.sub(r"\\item\[([^]]*)\]", r"\1:", s)
        s = re.sub(r"\\item\s*", "- ", s)
        s = re.sub(r"\\texttt\{([^}]*)\}", r"\1", s)
        s = s.replace(r"\[", "").replace(r"\]", "")
        s = s.replace(r"\mathrm", "")
        s = s.replace(r"\left", "").replace(r"\right", "")
        s = s.replace(r"\,", " ")
        s = s.replace("\\", "")
        s = s.replace("{", "").replace("}", "")
        lines.append(s if s.strip() else "")
    return lines


def wrap_lines(lines: list[str]) -> list[str]:
    wrapped = []
    for line in lines:
        if line.startswith("    "):
            wrapped.append(line)
        elif line.strip() == "":
            wrapped.append("")
        else:
            wrapped.extend(textwrap.wrap(line, width=92) or [""])
    return wrapped


def render_pdf(lines: list[str]) -> None:
    per_page = 48
    with PdfPages(OUT) as pdf:
        for start in range(0, len(lines), per_page):
            page_no = start // per_page + 1
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(
                0.06,
                0.95,
                "Atmosphere, Surface Water, and Hydrothermal Environment Architecture",
                fontsize=11,
                weight="bold",
                family="monospace",
            )
            y = 0.91
            for line in lines[start:start + per_page]:
                fig.text(0.06, y, line, fontsize=8.5, family="monospace", va="top")
                y -= 0.018
            fig.text(0.5, 0.03, f"Page {page_no}", fontsize=8, ha="center")
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    render_pdf(wrap_lines(latex_to_lines(SRC.read_text())))
    print(OUT)
