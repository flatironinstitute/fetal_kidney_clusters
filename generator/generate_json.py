
"""
Generate JSON data for kidney visualization from source data.
"""

import category3d

json_path = "../docs/data/all_categories.json"

def generate(
        data_path="./source/", 
        special_categories={9, 10, 11, 13}, 
        json_path=json_path,
        strides=4):
    print("Parsing", data_path)
    kidney = category3d.read_kidney_slices(path=data_path)
    print("Generating", json_path, "compressing", strides, "in each dimension.")
    category3d.all_categories_json(
        kidney, 
        to_json_path=json_path, 
        special_categories=special_categories,
        strides=strides)
    print("Open ../docs/index.html to view the generated visualization.")

if __name__ == "__main__":
    print(__doc__)
    generate()