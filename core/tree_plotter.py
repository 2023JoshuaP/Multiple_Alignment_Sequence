import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def _get_leaves(node):
    if node.is_leaf():
        return [node]
    return _get_leaves(node.left) + _get_leaves(node.right)

def plot_phylogenetic_tree(node, filepath, title="Árbol Filogenético", is_nj=False):
    """
    Dibuja un árbol (UPGMA o NJ) y lo guarda como imagen PNG.
    """
    plt.figure(figsize=(8, 6))
    
    leaves = _get_leaves(node)
    y_coords = {id(leaf): i for i, leaf in enumerate(leaves)}
    x_coords = {}
    
    if not is_nj:
        # UPGMA: node.height es la distancia desde las hojas.
        def get_max_height(n):
            if n.is_leaf(): return 0
            return max(n.height, get_max_height(n.left), get_max_height(n.right))
            
        max_h = get_max_height(node)
        
        def assign_coords_upgma(n):
            x_coords[id(n)] = max_h - n.height
            if n.is_leaf():
                return x_coords[id(n)], y_coords[id(n)]
            
            x_left, y_left = assign_coords_upgma(n.left)
            x_right, y_right = assign_coords_upgma(n.right)
            
            y_coords[id(n)] = (y_left + y_right) / 2.0
            return x_coords[id(n)], y_coords[id(n)]
            
        assign_coords_upgma(node)
    else:
        # NJ: node.height es la longitud de rama (branch length) hacia su padre.
        def assign_coords_nj_bottom_up(n):
            if n.is_leaf():
                return y_coords[id(n)]
            y_left = assign_coords_nj_bottom_up(n.left)
            y_right = assign_coords_nj_bottom_up(n.right)
            y_coords[id(n)] = (y_left + y_right) / 2.0
            return y_coords[id(n)]
            
        assign_coords_nj_bottom_up(node)
        
        x_coords[id(node)] = 0.0
        def assign_x_nj(n):
            if not n.is_leaf():
                x_coords[id(n.left)] = x_coords[id(n)] + n.left.height
                x_coords[id(n.right)] = x_coords[id(n)] + n.right.height
                assign_x_nj(n.left)
                assign_x_nj(n.right)
                
        assign_x_nj(node)
        
    def draw_lines(n):
        if not n.is_leaf():
            # Línea vertical
            plt.plot([x_coords[id(n)], x_coords[id(n)]], [y_coords[id(n.left)], y_coords[id(n.right)]], color='#2c3e50', lw=2.5, zorder=2)
            
            # Líneas horizontales a los hijos
            plt.plot([x_coords[id(n)], x_coords[id(n.left)]], [y_coords[id(n.left)], y_coords[id(n.left)]], color='#2c3e50', lw=2.5, zorder=2)
            plt.plot([x_coords[id(n)], x_coords[id(n.right)]], [y_coords[id(n.right)], y_coords[id(n.right)]], color='#2c3e50', lw=2.5, zorder=2)
            
            draw_lines(n.left)
            draw_lines(n.right)
            
    draw_lines(node)
    
    # Dibujar nodos como puntos
    all_x = list(x_coords.values())
    all_y = [y_coords[k] for k in x_coords.keys()]
    plt.scatter(all_x, all_y, color='#e74c3c', s=60, zorder=3, edgecolors='white', linewidths=1.5)
    
    # Etiquetas de hojas
    xs = list(x_coords.values())
    max_x = max(xs) if xs else 1
    min_x = min(xs) if xs else 0
    margin = (max_x - min_x) * 0.02 if max_x != min_x else 0.1
    
    for leaf in leaves:
        plt.text(x_coords[id(leaf)] + margin, y_coords[id(leaf)], leaf.label, va='center', ha='left', fontsize=12, fontweight='bold', color='#2c3e50')
        
    plt.title(title, fontsize=16, fontweight='bold', pad=20, color='#34495e')
    plt.xlabel("Distancia Evolutiva", fontsize=12, color='#7f8c8d')
    plt.yticks([])
    
    # Grid y estilos de ejes
    plt.grid(axis='x', linestyle='--', alpha=0.5, color='#bdc3c7', zorder=0)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_color('#7f8c8d')
    
    plt.xlim(min_x - margin*2, max_x + margin*15)
    
    plt.tight_layout()
    plt.savefig(filepath, format='png', dpi=200, bbox_inches='tight')
    plt.close()
    
    return filepath
