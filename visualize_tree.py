from tree import Tree, parse_test_data


def print_ascii(
    T, show_branch_lengths=True, show_internal_ids=False, sort_children=True
):
    root = T.find_root()

    def label(n, is_child=False):
        show_name = (
            (not n.children)
            or show_internal_ids
            or (not is_child and show_internal_ids)
        )
        base = str(n.name) if show_name else ""
        return base

    def rec(n, prefix="", is_last=True, from_parent=False):
        head = label(n, is_child=from_parent)
        if prefix == "":  # root line
            if head:
                print(head)
        else:
            branch = "└── " if is_last else "├── "
            bl = ""
            if from_parent and show_branch_lengths and n.branch_length is not None:
                bl = f" [{n.branch_length:.6g}]"
            name_part = head if head else ""
            # If no internal label requested, still show branch length after connector
            if not name_part and bl:
                print(prefix + branch + bl.strip())
            else:
                print(prefix + branch + (name_part + bl).strip())

        chs = n.children[:]
        if sort_children:
            chs.sort(key=lambda x: str(x.name))
        for i, ch in enumerate(chs):
            last = i == len(chs) - 1
            extender = "    " if is_last else "│   "
            rec(ch, prefix + extender, last, from_parent=True)

    rec(root)


if __name__ == "__main__":
    T = parse_test_data()
    print_ascii(T, show_branch_lengths=True, show_internal_ids=False)
