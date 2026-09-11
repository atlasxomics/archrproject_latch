"""Exercise real H5AD concatenation without importing plotting dependencies."""
import ast
import glob
import os
import logging
import tempfile
import unittest
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

namespace = dict(anndata=anndata, glob=glob, Path=Path, logging=logging)
module = ast.parse(Path('wf/features.py').read_text())
functions = [node for node in module.body if isinstance(node, ast.FunctionDef)
             and node.name in {'_combine_h5ad_files', '_ensure_anndata_root_encoding'}]
exec(compile(ast.Module(body=functions, type_ignores=[]), 'wf/features.py', 'exec'), namespace)


class DiskConcatTest(unittest.TestCase):
    def test_matches_in_memory_concat(self):
        for use_sparse in (False, True):
            with self.subTest(sparse=use_sparse), tempfile.TemporaryDirectory() as directory:
                previous = os.getcwd()
                os.chdir(directory)
                try:
                    objects = []
                    for i, genes in enumerate((['g2', 'g1', 'g3'], ['g3', 'g2', 'g4'])):
                        x = np.arange(6, dtype=np.float32).reshape(2, 3) + i
                        obj = anndata.AnnData(
                            sparse.csr_matrix(x) if use_sparse else x,
                            obs=pd.DataFrame({'sample': [str(i)] * 2}, index=[f'{i}a', f'{i}b']),
                            var=pd.DataFrame(index=genes),
                        )
                        obj.obsm['spatial'] = np.ones((2, 2)) * i
                        obj.layers['scores'] = x.copy()
                        obj.write_h5ad(f'{i}_g_converted.h5ad')
                        objects.append(obj)
                    # Exercise compatibility with older SeuratDisk root metadata.
                    with h5py.File('0_g_converted.h5ad', 'r+') as handle:
                        del handle.attrs['encoding-type']
                        del handle.attrs['encoding-version']
                    for file in glob.glob('*g_converted.h5ad'):
                        with h5py.File(file, 'r+') as handle:
                            for key in ('obsm/spatial', 'layers/scores'):
                                del handle[key].attrs['encoding-type']
                                del handle[key].attrs['encoding-version']
                    expected = anndata.concat([
                        anndata.read_h5ad(file)
                        for file in glob.glob("*g_converted.h5ad")
                    ])
                    actual = namespace['_combine_h5ad_files']('*g_converted.h5ad')
                    dense = lambda x: x.toarray() if sparse.issparse(x) else x
                    np.testing.assert_array_equal(dense(actual.X), dense(expected.X))
                    np.testing.assert_array_equal(actual.layers['scores'], expected.layers['scores'])
                    pd.testing.assert_frame_equal(actual.obs, expected.obs)
                    pd.testing.assert_frame_equal(actual.var, expected.var)
                    np.testing.assert_array_equal(actual.obsm['spatial'], expected.obsm['spatial'])
                    self.assertFalse(list(Path('.').glob('anndata_concat_*')))
                    with self.assertRaises(FileNotFoundError):
                        namespace['_combine_h5ad_files']('*missing.h5ad')
                finally:
                    os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
