import os
import argparse
from rawls.rawls import Rawls
import numpy as np

current_extension = '.rawls'
save_extension = '.png'

normals_file = 'normals'
zbuffer_file = 'zbuffer'

def extract(scene_folder, output):

    files = os.listdir(scene_folder)

    if scene_folder[-1] == '/':
        scene_folder = scene_folder[:-1]

    print(os.path.split(scene_folder))
    # build output data folder
    output_normals_scene_folder = os.path.join(output, normals_file, os.path.split(scene_folder)[-1])

    if not os.path.exists(os.path.join(output_normals_scene_folder)):
        os.makedirs(output_normals_scene_folder)

    output_zbuffer_scene_folder = os.path.join(output, zbuffer_file, os.path.split(scene_folder)[-1])

    if not os.path.exists(os.path.join(output_zbuffer_scene_folder)):
        os.makedirs(output_zbuffer_scene_folder)

    for filename in files:

        if normals_file in filename and current_extension in filename:

            file_path = os.path.join(scene_folder, filename)

            # build expected output files
            outfile_path = os.path.join(output_normals_scene_folder, filename.replace(current_extension, save_extension))
            
            print('extracts normals from', file_path, 'and save it into', outfile_path)

            # open file and normalize data
            normals_rawls = Rawls.load(file_path)
            # not necessary already normalized
            # normals_rawls.normalize() 

            normals_rawls.data = normals_rawls.data * 255
            
            # do not apply gamma conversion
            normals_rawls.save(outfile_path, gamma_convert=False)

        if zbuffer_file in filename and current_extension in filename:
            
            file_path = os.path.join(scene_folder, filename)

            # build expected output files
            outfile_path = os.path.join(output_zbuffer_scene_folder, filename.replace(current_extension, save_extension))
            
            print('extracts zbuffer from', file_path, 'and save it into', outfile_path)

            # open file and normalize data
            zbuffer_rawls = Rawls.load(file_path)

            zbuffer_rawls = zbuffer_rawls.normalize()
            zbuffer_rawls.data = zbuffer_rawls.data * 255

            # do not apply gamma conversion
            zbuffer_rawls.save(outfile_path, gamma_convert=False)

def main():

    parser = argparse.ArgumentParser(description="extract and save normals and zbuffer")

    parser.add_argument('--folder', type=str, help="folder which contains all scenes")
    parser.add_argument('--output', type=str, help="output data folder")

    args = parser.parse_args()

    p_folder = args.folder
    p_output = args.output

    for folder in os.listdir(p_folder):
        folder_path = os.path.join(p_folder, folder)
        extract(folder_path, p_output)

if __name__ == "__main__":
    main()
