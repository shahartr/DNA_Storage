
class TrieNode:
    def __init__(self):
        self.children = {}
        self.is_end_of_primer = False


class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, primer):
        node = self.root
        for char in primer:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.is_end_of_primer = True

    def count_tree_nodes(self):
        return self.count_tree_nodes_recursive(self.root)

    def count_tree_nodes_recursive(self, node):
        count = 1
        for child_node in node.children.values():
            count += self.count_tree_nodes_recursive(child_node)
        return count

    def count_leaves(self):
        return self.count_leaves_recursive(self.root)

    def count_leaves_recursive(self, node):
        if not node.children:
            return 1
        count = 0
        for child_node in node.children.values():
            count += self.count_leaves_recursive(child_node)
        return count

    # hamming distance is the number of positions at which the corresponding symbols are different
    #in our case we need every primer in the final set to have hamming distance of at least 6
    # basically we are searching (in tree that built from the primer added to the final set)
    # for a primer with hamming distance less than min_hamming
    # if we find such primer then we return False
    # (means that in the final 3set primer with hamming distance less than 6 is already present in the tree
    # and that why we cant add him to the final set
    def search_with_hamming_distance(self, node, primer, lvl, hamming_distance, min_hamming, histogram, found_invalid):
        # if hamming distance already more than 6 (6 mismatches found)) then valid
        if hamming_distance > min_hamming:
            if not found_invalid[0]:
                histogram[lvl-1] += 1  # increment at the level where it became valid
                found_invalid[0] = True
                return False

        # if end of primer check if the primer have hamming distance less than Min_Ham
        if lvl == len(primer):
            return node.is_end_of_primer and hamming_distance < min_hamming

        char = primer[lvl]
        for child_char, child_node in node.children.items():
            new_hamming_distance = hamming_distance + (1 if char != child_char else 0)

            if self.search_with_hamming_distance(child_node, primer, lvl + 1, new_hamming_distance, min_hamming, histogram, found_invalid):
                return True  # if found primer with hamming distance less than max_mismatches

        # if we've tried all children and haven't found a match
        if not found_invalid[0] and lvl == len(primer) - 1:
            histogram[lvl] += 1
            found_invalid[0] = True
        return False

    def is_valid_primer(self, primer, max_mismatches, histogram):
        found_invalid = [False]
        return not self.search_with_hamming_distance(self.root, primer, 0, 0, max_mismatches, histogram, found_invalid)
