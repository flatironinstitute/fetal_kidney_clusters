# fetal_kidney_clusters

<a href="https://flatironinstitute.github.io/fetal_kidney_clusters/">
A basic 3d visualization of fetal kidney cell type clusters.  View the visualization
at https://flatironinstitute.github.io/fetal_kidney_clusters/index.html.
</a>

To regenerate the visualization with new data clone the repository,
install new data in the `generatory/source` folder, and then run the
generator script:

```bash
% cd generator
% python3 generate_json.py 
```

The generator script requires `numpy` and `scipy.signal`.

