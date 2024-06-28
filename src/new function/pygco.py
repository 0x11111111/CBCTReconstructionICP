# encoding: utf-8
# module pygco
# from C:\Users\DYPro\.conda\envs\mesh_seg_net\lib\site-packages\pygco.cp39-win_amd64.pyd
# by generator 1.147
# no doc

# imports
import builtins as __builtins__ # <module 'builtins' (built-in)>
import numpy as np # C:\Users\DYPro\AppData\Roaming\Python\Python39\site-packages\numpy\__init__.py

# functions

def cut_from_graph(*args, **kwargs): # real signature unknown
    """
    Apply multi-label graphcuts to arbitrary graph given by `edges`.
    
        Parameters
        ----------
        edges: ndarray, int32, shape(n_edges, 2 or 3)
            Rows correspond to edges in graph, given as vertex indices.
            if edges is n_edges x 3 then third parameter is used as edge weight
        unary_cost: ndarray, int32, shape=(n_vertices, n_labels)
            Unary potentials
        pairwise_cost: ndarray, int32, shape=(n_labels, n_labels)
            Pairwise potentials for label compatibility
        n_iter: int, (default=5)
            Number of iterations
        algorithm: string, `expansion` or `swap`, default=expansion
            Whether to perform alpha-expansion or alpha-beta-swaps.
    """
    pass

def cut_simple(*args, **kwargs): # real signature unknown
    """
    Apply multi-label graphcuts to grid graph.
    
        Parameters
        ----------
        unary_cost: ndarray, int32, shape=(width, height, n_labels)
            Unary potentials
        pairwise_cost: ndarray, int32, shape=(n_labels, n_labels)
            Pairwise potentials for label compatibility
        n_iter: int, (default=5)
            Number of iterations
        algorithm: string, `expansion` or `swap`, default=expansion
            Whether to perform alpha-expansion or alpha-beta-swaps.
    """
    pass

def cut_simple_vh(*args, **kwargs): # real signature unknown
    """
    Apply multi-label graphcuts to grid graph.
    
        Parameters
        ----------
        unary_cost: ndarray, int32, shape=(width, height, n_labels)
            Unary potentials
        pairwise_cost: ndarray, int32, shape=(n_labels, n_labels)
            Pairwise potentials for label compatibility
        costV: ndarray, int32, shape=(width, height)
            Vertical edge weights
        costH: ndarray, int32, shape=(width, height)
            Horizontal edge weights
        n_iter: int, (default=5)
            Number of iterations
        algorithm: string, `expansion` or `swap`, default=expansion
            Whether to perform alpha-expansion or alpha-beta-swaps.
    """
    pass

# no classes
# variables with complex values

__loader__ = None # (!) real value is '<_frozen_importlib_external.ExtensionFileLoader object at 0x000002121D9A5700>'

__spec__ = None # (!) real value is "ModuleSpec(name='pygco', loader=<_frozen_importlib_external.ExtensionFileLoader object at 0x000002121D9A5700>, origin='D:\\\\Anaconda3\\\\envs\\\\MeshSegNet\\\\lib\\\\site-packages\\\\pygco.cp39-win_amd64.pyd')"

__test__ = {}

