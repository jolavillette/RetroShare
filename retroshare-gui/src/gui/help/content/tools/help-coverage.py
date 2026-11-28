#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Checks the manual of the GUI against the sources of the GUI.

Structure (an error makes the script fail):
  - every topic of en/contents.xml has its page, every page is a topic and is
    listed in content.qrc;
  - every link between pages points to an existing page;
  - every help topic given to a page of the GUI exists in the contents.

Coverage (a report; --strict makes it fail too):
  - the labels of the menu entries created in the C++ sources and of the
    buttons, check boxes, group boxes, tabs and actions of the .ui forms are
    looked for in the text of the manual. A label the manual never quotes is
    listed with the file that defines it.

Lines of tools/help-coverage.ignore are either a path prefix (relative to
retroshare-gui/src or to the repository root for the plugins), which skips
the files below it, or "label:" followed by a label to accept as missing.

Usage: help-coverage.py [--strict] [--quiet]
"""

import os
import re
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
CONTENT = os.path.dirname(HERE)
PAGES = os.path.join(CONTENT, "en")
GUI_SRC = os.path.normpath(os.path.join(CONTENT, "..", "..", ".."))
REPO = os.path.normpath(os.path.join(GUI_SRC, "..", ".."))

SOURCE_DIRS = [
    (os.path.join(GUI_SRC, "gui"), GUI_SRC),
    (os.path.join(REPO, "plugins", "VOIP"), REPO),
    (os.path.join(REPO, "plugins", "FeedReader"), REPO),
]

UI_TEXT_CLASSES = {
    "QPushButton", "QToolButton", "QCheckBox", "QRadioButton", "QCommandLinkButton",
}
UI_TITLE_CLASSES = {"QGroupBox"}

MENU_LABEL = re.compile(
    r'(?:addAction|QAction|addMenu)\s*\((?:[^;{}]{0,200}?)\btr\(\s*"((?:[^"\\]|\\.)*)"',
    re.S)
TOPIC_IN_CODE = re.compile(r'addPage\([^;]*,\s*"([a-z.]+)"\s*\)\s*;')
PAGE_LINK = re.compile(r'\]\(([^)#\s]+\.md)(?:#[^)]*)?\)')


def read(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        return handle.read()


def normalise(text):
    """Lower case, no accelerator, no markup, single spaces."""
    text = text.replace("&amp;", "&").replace("&&", "\0").replace("&", "").replace("\0", "&")
    text = re.sub(r"<[^>]+>", " ", text)
    text = text.replace("\\n", " ").replace('\\"', '"')
    text = re.sub(r"[*`_]", "", text)
    text = re.sub(r"\s+", " ", text).strip().lower()
    return text.rstrip(".:… ").strip()


def load_ignore():
    prefixes, labels = [], set()
    path = os.path.join(HERE, "help-coverage.ignore")
    if os.path.exists(path):
        for line in read(path).splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("label:"):
                labels.add(normalise(line[len("label:"):]))
            else:
                prefixes.append(line)
    return prefixes, labels


def check_structure():
    errors = []

    tree = ET.parse(os.path.join(PAGES, "contents.xml"))
    root = tree.getroot()
    topics = set()
    files = [root.get("html")]

    def walk(element, prefix):
        for child in element.findall("Topic"):
            topic = prefix + child.get("id")
            topics.add(topic)
            files.append(child.get("html"))
            walk(child, topic + ".")

    walk(root, "")

    on_disk = sorted(name for name in os.listdir(PAGES) if name.endswith(".md"))
    qrc = read(os.path.join(CONTENT, "content.qrc"))

    for name in sorted(set(files)):
        if name not in on_disk:
            errors.append("contents.xml: page %s does not exist" % name)
    for name in on_disk:
        if name not in files:
            errors.append("%s is not in contents.xml" % name)
        if "<file>en/%s</file>" % name not in qrc:
            errors.append("%s is not in content.qrc" % name)
        for target in PAGE_LINK.findall(read(os.path.join(PAGES, name))):
            if target not in on_disk:
                errors.append("%s: link to missing page %s" % (name, target))

    for source in ("gui/MainWindow.cpp", "gui/settings/rsettingswin.cpp"):
        for topic in TOPIC_IN_CODE.findall(read(os.path.join(GUI_SRC, source))):
            if topic not in topics:
                errors.append("%s: help topic %s is not in contents.xml" % (source, topic))

    return errors, on_disk


def ui_labels(path):
    labels = []
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return labels
    for widget in root.iter("widget"):
        cls = widget.get("class")
        for prop in widget.findall("property"):
            name = prop.get("name")
            if (cls in UI_TEXT_CLASSES and name == "text") or \
               (cls in UI_TITLE_CLASSES and name == "title"):
                string = prop.find("string")
                if string is not None and string.text:
                    labels.append(string.text)
        for attribute in widget.findall("attribute"):
            if attribute.get("name") == "title":
                string = attribute.find("string")
                if string is not None and string.text:
                    labels.append(string.text)
    for action in root.iter("action"):
        for prop in action.findall("property"):
            if prop.get("name") == "text":
                string = prop.find("string")
                if string is not None and string.text:
                    labels.append(string.text)
    return labels


def collect_labels(ignored_prefixes):
    found = {}
    for directory, base in SOURCE_DIRS:
        for folder, _, names in os.walk(directory):
            for name in sorted(names):
                path = os.path.join(folder, name)
                relative = os.path.relpath(path, base)
                if any(relative.startswith(prefix) for prefix in ignored_prefixes):
                    continue
                if name.endswith(".cpp"):
                    labels = MENU_LABEL.findall(read(path))
                elif name.endswith(".ui"):
                    labels = ui_labels(path)
                else:
                    continue
                for label in labels:
                    key = normalise(label)
                    # A label without two letters in a row is an icon button.
                    if not re.search(r"[a-z]{2}", key):
                        continue
                    found.setdefault(relative, [])
                    if (key, label) not in found[relative]:
                        found[relative].append((key, label))
    return found


def main():
    strict = "--strict" in sys.argv
    quiet = "--quiet" in sys.argv

    errors, pages = check_structure()
    for error in errors:
        print("ERROR  " + error)

    manual = normalise(" ".join(read(os.path.join(PAGES, name)) for name in pages))
    ignored_prefixes, ignored_labels = load_ignore()
    labels = collect_labels(ignored_prefixes)

    total = missing = 0
    for relative in sorted(labels):
        absent = []
        for key, label in labels[relative]:
            total += 1
            if key in manual or key in ignored_labels:
                continue
            missing += 1
            absent.append(label)
        if absent and not quiet:
            print(relative)
            for label in absent:
                print("    " + label.replace("\n", " "))

    print("%d pages, %d structure errors, %d labels checked, %d not in the manual"
          % (len(pages), len(errors), total, missing))

    if errors or (strict and missing):
        sys.exit(1)


if __name__ == "__main__":
    main()
