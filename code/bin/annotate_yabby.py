#!/usr/bin/env python

from ete3 import Tree
import os 
import re 
# handle arguments
import argparse


parser = argparse.ArgumentParser(description="automatic coloration and keep label in figtree for yabby trees")

parser.add_argument("tree", help="path towards tree", type=str)

args = parser.parse_args()
tree_path = args.tree

print(tree_path)

path_name = os.path.dirname(os.path.abspath(tree_path))

print(path_name)

rename_tree = Tree(tree_path, format=1)

col_dict= {"YABA":"#ffe100","YABB":"#00b2ff","YABC":"#00ff0f","YABD":"#2300ff","FIL":"#c800ff","CRC":"#ff8500","INO":"#ff00d2","YAB5":"#ff001c","YAB2":"#0aa913"}


for leaf in rename_tree.traverse():
    for yab in col_dict.keys():
        if yab in leaf.name.split("_"):
            leaf.add_features(yabby=yab) ## add a yabby feature to the leaves
            leaf.add_features(color=col_dict[yab])
for node in rename_tree.traverse():
    result = re.search(r'\["(\d+)"\]', node.name) ##look for bootstrap in between quotes
    if not result is None:
        node.add_features(label = result.group(1)) ##add bootstrap as an attribute

for i in col_dict:
    for node in rename_tree.get_monophyletic(values=[i], target_attr="yabby"): #get all the yabby subtrees
        for descendent in node.get_descendants():
            descendent.add_features(color = col_dict[i]) #annotate all the nodes in a subtree with the yab name
            node.add_features(color = col_dict[i]) 


for node in rename_tree.traverse():
    if ("label" in node.features) and ("color" in node.features):
        comment = node.name.split("[")[0]+"[&label="+str(node.label)+",!color="+node.color+"]"
    elif ("color" in node.features) and ("label" not in node.features):
        comment = node.name+"[&!color="+node.color+"]"
    elif ("color" not in node.features) and ("label" in node.features):
        comment = node.name.split("[")[0]+"[&label="+node.label+"]"
    else: 
        comment = node.name+"[&!color=#000000]"
    node.add_features(annotation = comment)

with open(path_name+ "/internal_nodes_colors.txt", "w") as wf:
    for node in rename_tree.traverse():
        wf.write(node.name.split("[")[0]+"\t"+node.annotation+"\n")