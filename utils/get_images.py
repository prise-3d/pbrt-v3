import os
import argparse
from rawls.rawls import Rawls

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
    output_scene_folder = os.path.join(output, os.path.split(scene_folder)[-1])

    if not os.path.exists(os.path.join(output_scene_folder)):
        os.makedirs(output_scene_folder)

    for filename in files:

        if (normals_file not in filename and zbuffer_file not in filename) and current_extension in filename:

            file_path = os.path.join(scene_folder, filename)

            # build expected output files
            outfile_path = os.path.join(output_scene_folder, filename.replace(current_extension, save_extension))
            
            print('extracts image from', file_path, 'and save it into', outfile_path)

            # open file and normalize data
            images_rawls = Rawls.load(file_path)

            # do not apply gamma conversion
            images_rawls.save(outfile_path)

def main():

    parser = argparse.ArgumentParser(description="extract and save normals and zbuffer")

    parser.add_argument('--folder', type=str, help="folder which contains normals and zbuffer data (rawls format)")
    parser.add_argument('--output', type=str, help="output data folder")

    args = parser.parse_args()

    p_folder = args.folder
    p_output = args.output

    extract(p_folder, p_output)

if __name__ == "__main__":
    main()
