# class to create a tree structure representing
class Node:
    """
    label: name of the node

    """
    def __init__(self, label):
        self.label = label
        self.distance = 0
        self.children = {}
        self.sequence = ''

class NewickTree:
    def __init__(self, label):