import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input')
parser.add_argument('-p', '--path')
parser.add_argument('-v', '--value')
parser.add_argument('-o', '--output')
args = parser.parse_args()

# Assume XML conformant parameters
try:
    root = ET.parse(args.input).getroot()
# Old non-conformant parameters
except ET.ParseError:
    with open(args.input) as f:
        next(f) # 
        # introduce missing XML root
        xml_string = '<quicc>\n' + f.read() + '</quicc>'

    root = ET.fromstring(xml_string)

# Find requested node and update value
node = root.find(args.path)
node.text = args.value

# Create new XML file
tree = ET.ElementTree(root)
tree.write(args.output, encoding='utf-8', xml_declaration=True)
