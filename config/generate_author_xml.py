import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
# import codecs # pip install latexcodec
# import latexcodec
from datetime import datetime
from pylatexenc.latex2text import LatexNodes2Text # pip install pylatexenc

# --- INPUT INFO: Adjust as needed --- #
# Paths to your LaTeX input files
AUTHORS_TEX = "authors.tex"
INSTITUTIONS_TEX = "institutions.tex"
OUTPUT_XML = "author.xml"

# Set your arXiv publication reference
publication_url = "http://arXiv.org/abs/2509.05758"

# --- Set date and time --- #
# Current timestamp in the format: YYYY-MM-DD_HH:MM
creation_date = datetime.now().strftime("%Y-%m-%d_%H:%M")

# --- Set latex decoding function --- #
# def decode_latex(text: str) -> str:
#     if not text:
#         return ""
    # text = text.replace(r'`{', '\\`{')

#     # Step 1: Ensure input is a string and encode to bytes
#     try:
#         bytes_text = text.encode('utf-8')
#     except AttributeError:
#         # Already bytes
#         bytes_text = text

#     # Step 2: Decode using latexcodec
#     try:
#         decoded = bytes_text.decode('latex')
#     except Exception as e:
#         print(f"Error decoding LaTeX: {e}")
#         return text  # Return raw if decoding fails

#     # Step 3: Remove braces around single characters (e.g., {é} → é)
#     cleaned = re.sub(r'\{([^\{\}])\}', r'\1', decoded)

#     # Optional: Collapse multiple spaces and strip
#     return cleaned.strip()

def decode_latex(text: str) -> str:
    if not text:
        return ""

    text = text.replace(r'`{', '\\`{')

    # Correct usage: create a persistent LatexNodes2Text instance
    latex_converter = LatexNodes2Text()
    decoded = latex_converter.latex_to_text(text)

    # Optional cleanup: remove stray braces (e.g., around single letters)
    cleaned = re.sub(r'\{([^\{\}])\}', r'\1', decoded)

    return cleaned.strip()

# --- Parse institutions.tex ---

with open(INSTITUTIONS_TEX, "r", encoding="utf-8") as f:
    tex = f.read()

# Matches: \newcommand*{\DUKE}{Duke University...}
aff_pattern = re.compile(r'\\newcommand\*?{\\(\w+)}\s*{((?:[^{}]|\\{|\\})+)}')#r'\\newcommand\*?\s*{\\(\w+)}\s*{([^}]+)}')
aff_pattern = re.compile(r'\\newcommand\*?{\\(\w+)}\s*{((?:[^{}\\]|\\.|{[^{}]*})+)}')
# Matches: \newcommand*{\DUKEindex}{6}
index_pattern = re.compile(r'\\newcommand\*?\s*{\\(\w+)index}\s*{(\d+)}')

# Build macro -> index map
indices = {match.group(1): match.group(2) for match in index_pattern.finditer(tex)}

# Build macro -> {name, id}
aff_map = {}
for match in aff_pattern.finditer(tex):
    macro, name = match.groups()
    if macro.endswith("index"):
        continue  # skip DUKEindex etc.
    name = decode_latex(name.strip())
    start_idx = 1 if name.startswith("INFN") else 0
    orgName = ", ".join(name.split(", ")[:start_idx+1])
    orgAddress = ", ".join(name.split(", ")[start_idx+1:])
    aff_map[macro] = {"name": orgName, "address":orgAddress, "id": indices.get(macro, "?")}

# --- Parse authors.tex ---

with open(AUTHORS_TEX, "r", encoding="utf-8") as f:
    lines = f.readlines()

authors = []
current_author = {}
for line in lines:
    line = line.strip()
    if line.startswith(r"\author"):
        if current_author:
            authors.append(current_author)
        name_match = re.search(r"{(.+?)}", line)
        name = name_match.group(1) if name_match else ""
        name = name.replace("~", " ")  # Replace tilde with space
        name = decode_latex(name)
        current_author = {
            "name": name,
            "email": "",
            "affiliations": [],
            "alt": None
        }
    elif line.startswith(r"\email"):
        email_match = re.search(r"{(.+?)}", line)
        current_author["email"] = email_match.group(1) if email_match else ""
    elif line.startswith(r"\affiliation"):
        aff_match = re.search(r"{\\(\w+)}", line)
        if aff_match:
            current_author["affiliations"].append(aff_match.group(1))
    elif line.startswith(r"\altaffiliation"):
        alt_match = re.search(r"{(.+?)}", line)
        if alt_match:
            alt_text = alt_match.group(1).replace(r'\NOW', '').replace("~", " ")
            current_author["alt"] = alt_text.strip()

if current_author:
    authors.append(current_author)

# --- Create XML --- #

root = ET.Element(
    "collaborationauthorlist",
    {
        "xmlns:foaf": "http://xmlns.com/foaf/0.1/",
        "xmlns:cal": "http://inspirehep.net/info/HepNames/tools/authors_xml/"
    }
)

ET.SubElement(root, "cal:creationDate").text = creation_date
ET.SubElement(root, "cal:publicationReference").text = publication_url

# Add collaboration elements
collaborations = ["CLAS"]
collaborationid = None #NOTE: Assumee only one collaboration for now
collaborations_el = ET.SubElement(root, "cal:collaborations")
for idx, cname in enumerate(collaborations):
    collaborationid_el = "c"+str(idx + 1)
    collaboration_el = ET.SubElement(collaborations_el, "cal:collaboration", id=collaborationid_el)
    ET.SubElement(collaboration_el, "foaf:name").text = cname
    collaborationid = collaborationid_el

# Add <author> elements
authors_el = ET.SubElement(root, "cal:authors")
for author in authors:
    author_el = ET.SubElement(authors_el, "foaf:Person")
    collaboration_el = ET.SubElement(author_el, "cal:authorCollaboration", collaborationid=collaborationid)
    name_el = ET.SubElement(author_el, "cal:authorNamePaper")
    gname_el = ET.SubElement(author_el, "cal:authorNamePaperGiven")
    fname_el = ET.SubElement(author_el, "cal:authorNamePaperFamily")
    name_el.text = author["name"]
    gname_el.text = author["name"].split(" ")[0]
    fname_el.text = " ".join(author["name"].split(" ")[1:])

    # if author["email"]:
    #     email_el = ET.SubElement(author_el, "email")
    #     email_el.text = author["email"]

    author_affs_el = ET.SubElement(author_el, "cal:authorAffiliations")
    for aff_macro in author["affiliations"]:
        aff_info = aff_map.get(aff_macro)
        if aff_info and aff_info["id"].isdigit():

            aff_el = ET.SubElement(author_affs_el, "cal:authorAffiliation", organizationid=aff_info["id"])
            # aff_el.text = aff_info["name"]
        else:
            print(f"⚠️ Warning: Unknown or unindexed affiliation macro: {aff_macro}")

    # if author["alt"]:
    #     note_el = ET.SubElement(author_el, "note")
    #     note_el.text = author["alt"]

# Add <institutions> section
institutions_el = ET.SubElement(root, "cal:organizations")
for macro, data in sorted(
    ((m, d) for m, d in aff_map.items() if d["id"].isdigit()),
    key=lambda x: int(x[1]["id"])
):
    inst_el = ET.SubElement(institutions_el, "foaf:Organization", id=data["id"])
    inst_name_el = ET.SubElement(inst_el, "foaf:name")
    inst_orgName_el = ET.SubElement(inst_el, "cal:orgName")
    inst_orgAddress_el = ET.SubElement(inst_el, "cal:orgAddress")
    inst_name_el.text = data["name"]
    inst_orgName_el.text = data["name"]
    inst_orgAddress_el.text = data["address"]

# --- Pretty print and write to file ---

def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_str = ET.tostring(elem, encoding="utf-8")
    reparsed = minidom.parseString(rough_str)
    return reparsed.toprettyxml(indent="  ")

with open(OUTPUT_XML, "w", encoding="utf-8") as f:
    pretty_xml = prettify(root)
    pretty_xml_no_decl = "\n".join(pretty_xml.split("\n")[1:])  # Remove XML declaration
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n') # Replace XML declaration
    f.write('<!DOCTYPE collaborationauthorlist SYSTEM "author.dtd">\n')
    f.write(pretty_xml_no_decl)

print(f"✅ Generated '{OUTPUT_XML}' with {len(authors)} authors and {len(aff_map)} affiliations.")
