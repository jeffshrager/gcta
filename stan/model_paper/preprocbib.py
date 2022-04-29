#!/bin/env python3
import json
# from pyzotero import zotero
from pybtex import database as db

# # download zotero keys
# kwargs = json.load(open("pyzotero_inputs.json"))
# zot = zotero.Zotero(**kwargs)

# read in data from zotero output
bd = db.parse_file("gcta.bib")

entries = []
# append setting for name-year format in natbib
entries.append(("label", db.Entry("settings", fields=[("options", "nameyear")])))
# add all entries with both names and years
for key, entry in bd.entries.items():
    if entry.type == "settings":
        entries.append((key, entry))
    if (entry.fields.get("year") is None) or (entry.persons.get("author") is None):
        continue
    entries.append((key, entry))

# output cleaned data
clean_bd = db.BibliographyData(entries=entries,
                               preamble=bd.preamble,
                               wanted_entries=bd.wanted_entries,
                               min_crossrefs=bd.min_crossrefs)
clean_bd.to_file("gcta_clean.bib")    

        
