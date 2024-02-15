class graph_obj:
    def __init__(self):
        self.nodes = []
        self.edges = []
        self.node_label_size = None
        self.node_label_color= None
        self.edge_label_size = None
        self.edge_label_color= None
    def set_metadata(self, node_label_size, node_label_color, edge_label_size, edge_label_color):
        self.node_label_size = node_label_size
        self.node_label_color= node_label_color
        self.edge_label_size = edge_label_size
        self.edge_label_color= edge_label_color