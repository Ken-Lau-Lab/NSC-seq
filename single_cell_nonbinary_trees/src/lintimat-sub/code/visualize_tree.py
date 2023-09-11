# read the input newick tree and create a visualization
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, faces
import os
import matplotlib.pyplot as plt
import copy


# read all meta-data
# cell --> number
# cell --> type --> color
# cell --> mutation --> names


# read the data


def get_mapping(
    cell_map_mutation_path,
    cell_map_type_path,
    cell_type_map_color_path,
    mutation_map_names_path,
):
    cell_map_mutation = pd.read_csv(cell_map_mutation_path, sep="\t", header=None)
    cell_map_mutation.append(
        [["normal", "-".join(["NONE"] * len(cell_map_mutation.iloc[0, 1].split("-")))]]
    )
    cell_map_mutation = {row[0]: row[1] for i, row in cell_map_mutation.iterrows()}
    cells = list(cell_map_mutation.keys())

    cell_map_type = pd.read_csv(cell_map_type_path, sep="\t", header=0, index_col=False)
    cell_map_type.set_index(pd.Index(cell_map_type["Cells"]), inplace=True)
    cell_map_type = cell_map_type.loc[cells]
    cell_map_type = {row[0]: row[1] for _, row in cell_map_type.iterrows()}

    # cell name and its type
    cell_type_map_color = pd.read_csv(
        cell_type_map_color_path, sep="\t", header=0, index_col=0
    )
    cell_type_map_color = {
        row[0]: "#" + row[1] for _, row in cell_type_map_color.iterrows()
    }

    mutation_map_names = pd.read_csv(mutation_map_names_path, header=None, sep="\t")[
        0
    ].to_list()

    return cell_map_mutation, cell_map_type, cell_type_map_color, mutation_map_names


def get_mutation_all(mutation, mutation_map_names):
    mutation_all = list()
    if mutation == "":
        return [set(["NONE"]) for i in range(len(mutation_map_names))]
    for mut in mutation.split("-"):
        if mut == "NONE":
            mutation_all.append(set(["NONE"]))
        else:
            mutation_all.append(set([mutation_map_names[int(mut)]]))
    return mutation_all


def add_attributes(
    tree, cell_map_mutation, cell_map_type, cell_type_map_color, mutation_map_names
):
    if not hasattr(add_attributes, "thisid"):
        add_attributes.thisid = 0  # it doesn't exist yet, so initialize it
    add_attributes.thisid += 1
    name = tree.name
    cell_type = cell_map_type.get(name, None)
    nodecolor = cell_type_map_color.get(cell_type, "#FFFFFF")
    mutation = cell_map_mutation.get(name, "")
    mutation_all = get_mutation_all(mutation, mutation_map_names)

    tree.add_features(
        **{
            "id": "i" + str(add_attributes.thisid),
            "cell_type": cell_type,
            "nodecolor": nodecolor,
            "mutation": mutation,
            "mutation_all": mutation_all,
            "mutation_diff": "",
        }
    )

    for child in tree.children:
        add_attributes(
            child,
            cell_map_mutation,
            cell_map_type,
            cell_type_map_color,
            mutation_map_names,
        )


# as the last branch should not have mutation adding dummpy parent to all leafs
def add_parent_to_leafs(tree):
    if tree.is_leaf():
        internal_node = Tree()
        internal_node.name = "added_node"
        internal_node.id = "-1"
        tree.up.add_child(internal_node)

        leaf = tree.detach()
        internal_node.add_child(leaf)
    else:
        children = tree.children.copy()
        for child in children:
            add_parent_to_leafs(child)


def build_mutations_set(tree):
    score = 0
    if tree is not None:
        if len(tree.children) != 0:
            for child in tree.children:
                score += build_mutations_set(child)
            if len(tree.children) == 1:
                tree.mutation_all = copy.deepcopy(tree.children[0].mutation_all)
            elif len(tree.children) == 2:
                for i in range(len(tree.mutation_all)):
                    if (
                        len(
                            tree.children[0]
                            .mutation_all[i]
                            .intersection(tree.children[1].mutation_all[i])
                        )
                        == 0
                    ):
                        tree.mutation_all[i] = (
                            tree.children[0]
                            .mutation_all[i]
                            .union(tree.children[1].mutation_all[i])
                        )
                        score += 1
                    else:
                        tree.mutation_all[i] = (
                            tree.children[0]
                            .mutation_all[i]
                            .intersection(tree.children[1].mutation_all[i])
                        )
            else:
                print("error recieved non exptected tree > 2 children  ")
    return score


def finalize_mutation(tree, is_root=True, parent_mutation=[]):
    if tree is not None:
        mutation_all = []
        if is_root:
            for mut in tree.mutation_all:
                mutation_all.append(["NONE"])
        else:
            for i in range(len(parent_mutation)):
                if len(parent_mutation[i]) == 0 and len(tree.mutation_all[i]) != 0:
                    print(
                        "------error----- parent has no mutation but child has",
                        tree.id,
                        tree.up.id,
                    )
                if parent_mutation[i][0] in tree.mutation_all[i]:
                    mutation_all.append(parent_mutation[i].copy())
                else:
                    mutation_all.append([tree.mutation_all[i].pop()])

        tree.mutation_all = mutation_all
        tree.add_features(
            **{"mutation_str": str([mut for mut in tree.mutation_all if len(mut) != 0])}
        )
        for child in tree.children:
            finalize_mutation(child, False, tree.mutation_all)


def order_tree(tree):
    if tree is not None:
        tree.children = sorted(
            tree.children,
            key=lambda x: [
                len(x.get_leaf_names()),
                len(x.get_descendants()),
                *[leaf.cell_type for leaf in x.get_leaves()],
            ],
        )
        for child in tree.children:
            order_tree(child)


def add_mutation_diff(tree):
    if tree is not None:
        for child in tree.children:
            diff = ""
            for i in range(len(child.mutation_all)):
                if child.mutation_all[i] != tree.mutation_all[i]:
                    if child.mutation_all[i] != ["NONE"]:
                        diff = diff + "\n" + str(child.mutation_all[i])
                    else:
                        diff = diff + "\n" + str(tree.mutation_all[i])
            if diff == "":
                diff = "NONE"
            child.add_feature("mutation_diff", diff)
            add_mutation_diff(child)


def cal_score(tree):
    score = 0
    if len(tree.children) != 0:
        par_mut = tree.mutation_all
        for child in tree.children:
            chi_mut = child.mutation_all
            for i in range(len(par_mut)):
                if par_mut[i] != chi_mut[i]:
                    score += 1
        for child in tree.children:
            score += cal_score(child)
    return score


def remove_non_mutation_branches(tree, is_first=True):
    # if is_first:
    #     print('removing ', end=' ')
    if tree.children != []:  # if not leaf
        check_tree_again = False
        for child in tree.children:
            if child.is_leaf():
                continue

            def same_mutation(parent_mut, child_mut):
                for i in range(len(parent_mut)):
                    if parent_mut[i] != child_mut[i]:
                        return False
                return True

            if same_mutation(tree.mutation_all, child.mutation_all):
                # print(child.id, end=',')
                if not child.is_leaf():  # internal node
                    child.delete()
                    check_tree_again = True
        if check_tree_again:
            if tree.is_root():
                remove_non_mutation_branches(tree, False)
            else:
                remove_non_mutation_branches(tree.up, False)
        else:
            for child in tree.children:
                remove_non_mutation_branches(child, False)
    if is_first:
        print()


# make the normal node as the root
def add_rootDist(tree, rootDist):
    if tree is None:
        return
    tree.rootDist = rootDist
    for child in tree.children:
        add_rootDist(child, rootDist + 1)


def make_normal_as_root(tree):
    normal = tree.search_nodes(name="normal")
    tree.set_outgroup(normal[0])
    for child in tree.children:
        try:
            print(child.id)
        except Exception:
            # print("got exception in normal as root child.id")
            for feat in tree.features:
                child.add_feature(feat, tree.__getattribute__(feat))
            child.id = "9999"


def same_level_leaf(tree):
    def equal_distance(tree):
        if tree.is_leaf():
            return
        for child in tree.children:
            equal_distance(child)
        sub_tree_depths = []
        for child in tree.children:
            farthest_leaf, farthest_leaf_distance = child.get_farthest_leaf()
            sub_tree_depths.append(farthest_leaf_distance)
        max_far = max(sub_tree_depths)
        for i, child in enumerate(tree.children):
            child.dist = max_far - sub_tree_depths[i] + 1

    def leaf_distance_1(tree):
        if tree.is_leaf():
            return
        leaf_dis = 1
        internal_node = Tree()
        internal_node.mutation_all = tree.mutation_all
        leafs_to_rearrange = []
        # get the nodes to remove
        for child in tree.children:
            if child.is_leaf() and child.dist > 1.0:
                leaf_dis = child.dist
                leafs_to_rearrange.append(child)
                child.dist = 1
        if leaf_dis != 1:
            new_children = []
            # remove the nodes
            for node in leafs_to_rearrange:
                new_children.append(node.detach())
            # add the nodes back to an internal node
            for node in new_children:
                internal_node.add_child(node)
            internal_node.dist = leaf_dis - 1
            tree.children[0:0] = [internal_node]
            internal_node.up = tree
            node_style = NodeStyle()
            node_style["hz_line_type"] = 2
            internal_node.set_style(node_style)
        for child in tree.children:
            leaf_distance_1(child)

    equal_distance(tree)
    # leaf_distance_1(tree)


def style(tree, cell_type_map_color, cwd, is_circular=True, is_same_level_leaf=True):
    color_annotation(cell_type_map_color, cwd)
    ts = TreeStyle()
    ts.allow_face_overlap = True

    if is_circular:
        ts.mode = "c"
    ts.show_leaf_name = False
    img_face = faces.ImgFace(cwd + "legend.png")  # , width = 100)
    ts.title.add_face(img_face, column=0)
    for node in tree.traverse():
        if node.is_leaf():
            nstyle = NodeStyle()
            nstyle["fgcolor"] = node.nodecolor
            nstyle["size"] = 20
            nstyle["shape"] = "circle"
            nstyle["hz_line_color"] = node.nodecolor
            nstyle["hz_line_width"] = 10
            node.set_style(nstyle)
        else:
            # pass
            cell_type_pie(node, cell_type_map_color, cwd)
        # mutation_diff = TextFace(node.mutation_diff, fgcolor='black',fsize=50)
        # node.add_face(mutation_diff, column=0, position="branch-top")
    if is_same_level_leaf:
        same_level_leaf(tree)
    return ts


def is_binary_tree(tree):
    if (
        len(tree.children) == 2
        and is_binary_tree(tree.children[0])
        and is_binary_tree(tree.children[1])
    ) or len(tree.children) == 0:
        return True
    else:
        return False


def max_depth(tree):
    if tree.children == []:
        return 0
    return 1 + max(max_depth(tree.children[0]), max_depth(tree.children[1]))


def print_attributes(tree):
    if tree is not None:
        print(tree.id)
        print("   {} {}".format("name", tree.name))
        for feat in tree.features:
            print("   {} {}".format(feat, tree.__getattribute__(feat)))

        for child in tree.children:
            print_attributes(child)


def is_one_to_zero(tree):
    if tree is not None:
        for child in tree.children:
            for i in range(len(tree.mutation_all)):
                if tree.mutation_all[i] != ["NONE"] and child.mutation_all == ["NONE"]:
                    print(
                        "--------------- one to zero ------------form id {} to id {}".format(
                            tree.id, child.id
                        )
                    )
                    raise Exception("There is one to zero conversion of mutation")
            is_one_to_zero(child)


def color_annotation(cell_type_map_color, cwd):
    import matplotlib.pyplot as plt

    # Create a plot with some data and legends
    fig, ax = plt.subplots()
    for key, val in cell_type_map_color.items():
        ax.plot([1], [1], label=key, color=val, linestyle="", marker="o")
    ax.legend()
    fig2, ax2 = plt.subplots()
    ax2.axis("off")
    legend = ax2.legend(*ax.get_legend_handles_labels(), loc="center", frameon=False)
    height = 1
    aspect_ratio = fig2.get_size_inches()[0] / fig2.get_size_inches()[1]
    width = height * aspect_ratio
    fig2.set_size_inches(width, height)
    fig2.savefig(cwd + "legend.png", dpi=500, transparent=True, bbox_inches="tight")
    print('legend saved to : ',cwd + "legend.png")
    plt.close(fig)
    plt.close(fig2)


def cell_type_pie(node, cell_type_map_color, cwd):
    if not hasattr(cell_type_pie, "pie_idx"):
        cell_type_pie.pie_idx = 0
    data = [leaf.cell_type for leaf in node.get_leaves()]
    if data[0] == "":
        print("node name is : ", node.name)
    data = pd.Series(data).value_counts()
    colors = [cell_type_map_color[val] for val in data.index]

    # Data
    labels = data.index
    sizes = data.values

    # Create pie chart
    fig, ax = plt.subplots()
    ax.pie(sizes, colors=colors, startangle=90)  # labels=labels,autopct='%1.1f%%',

    # Add count to legend
    # total = sum(sizes)
    # legend_labels = [f'({size}, {size/total*100:.1f}%)' for label, size in zip(labels, sizes)]
    # ax.legend(legend_labels, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    height = 1
    aspect_ratio = fig.get_size_inches()[0] / fig.get_size_inches()[1]
    width = height * aspect_ratio
    fig.set_size_inches(width, height)
    # plt.pie(data, colors=colors, autopct='%1.1f%%')
    plt.savefig(cwd + "temp_pie" + str(cell_type_pie.pie_idx) + ".png")  # dpi=500
    plt.close()
    img_face = faces.ImgFace(
        cwd + "temp_pie" + str(cell_type_pie.pie_idx) + ".png"
    )  # , width = 10)
    cell_type_pie.pie_idx += 1
    node.add_face(img_face, column=0, position="branch-top")


def del_pie_local(cwd):
    for pie_idx in range(cell_type_pie.pie_idx):
        os.remove(cwd + "temp_pie" + str(pie_idx) + ".png")


def convert_to_non_binary(
    tree,
    cell_map_mutation_path,
    cell_map_type_path,
    cell_type_map_color_path,
    mutation_map_names_path,
    cwd
):
    (
        cell_map_mutation,
        cell_map_type,
        cell_type_map_color,
        mutation_map_names,
    ) = get_mapping(
        cell_map_mutation_path,
        cell_map_type_path,
        cell_type_map_color_path,
        mutation_map_names_path,
    )
    add_parent_to_leafs(tree)
    add_attributes(
        tree, cell_map_mutation, cell_map_type, cell_type_map_color, mutation_map_names
    )
    print("parsimony_score ", build_mutations_set(tree))
    finalize_mutation(tree)
    # print_attributes(tree)
    remove_non_mutation_branches(tree)
    print("parsimony score cal_score:", cal_score(tree))
    make_normal_as_root(tree)
    add_mutation_diff(tree)
    order_tree(tree)
    is_one_to_zero(tree)
    # tree_sll = tree.copy(method="deepcopy")

    # ts_sll = visualize_tree.style(
    #     tree_sll,
    #     cell_type_map_color,
    #     cwd + "../output/",
    #     is_circular=True,
    #     is_same_level_leaf=True,
    # )
    # tree_sll.render(
    #     cwd + "../output/tree_ete3_sll.pdf", w=10000, units="px", dpi=500, tree_style=ts_sll
    # )

    tree_style = style(
        tree,
        cell_type_map_color,
        cwd + "../output/",
        is_circular=True,
        is_same_level_leaf=False,
    )
    return tree, tree_style
