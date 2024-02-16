class graph_obj:
    def __init__(self):
        self.graph_gjgf = {}
        self.graph_gjgf['graph'] = {}
        self.graph_gjgf['graph']['metadata'] = {}
        self.graph_gjgf['graph']['directed'] = True
        self.graph_gjgf['graph']['nodes'] = []
        self.graph_gjgf['graph']['edges'] = []
    def set_metadata(self, node_label_size, node_label_color, edge_label_size, edge_label_color):
        self.graph_gjgf['graph']['metadata']['node_label_size'] = node_label_size
        self.graph_gjgf['graph']['metadata']['node_label_color']= node_label_color
        self.graph_gjgf['graph']['metadata']['edge_label_size'] = edge_label_size
        self.graph_gjgf['graph']['metadata']['edge_label_color']= edge_label_color
    def add_node(self, id, label, color, size, shape, label_color, label_size, border_color, border_size):
        node = {}
        node['id'] = id
        node['label'] = label
        node['metadata'] = {}
        node['metadata']['color'] = color
        node['metadata']['size']  = size
        node['metadata']['shape'] = shape
        node['metadata']['opacity'] = 0.7
        node['metadata']['label_color'] = label_color
        node['metadata']['label_size'] = label_size
        node['metadata']['border_color'] = border_color
        node['metadata']['border_size'] = border_size
        self.graph_gjgf['graph']['nodes'].append(node)