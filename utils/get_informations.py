import os
import argparse
from rawls.rawls import Rawls

current_extension = '.rawls'
save_extension = '.png'

normals_file = 'normals'
zbuffer_file = 'zbuffer'

def extract(scene_folder, output):

    files = os.listdir(scene_folder)

    # build output data folder
    output_scene_folder = os.path.join(output, scene_folder)

    if not os.path.exists(os.path.join(output_scene_folder)):
        os.makedirs(output_scene_folder)

    for filename in files:

        if normals_file in filename and current_extension in filename:

            file_path = os.path.join(scene_folder, filename)

            # build expected output files
            outfile_path = os.path.join(output_scene_folder, filename.replace(current_extension, save_extension))
            
            # open file and normalize data
            normals_rawls = Rawls.load(file_path)
            normals_rawls.normalize()

            normals_rawls.data = normals_rawls.data * 255

            # do not apply gamma conversion
            normals_rawls.save(outfile_path, gamma_convert=False)

        if zbuffer_file in filename and current_extension in filename:
            
            file_path = os.path.join(scene_folder, filename)

            # build expected output files
            outfile_path = os.path.join(output_scene_folder, filename.replace(current_extension, save_extension))
            
            # open file and normalize data
            zbuffer_rawls = Rawls.load(file_path)
            zbuffer_rawls.normalize()

            zbuffer_rawls.data = zbuffer_rawls.data * 255

            # do not apply gamma conversion
            zbuffer_rawls.save(outfile_path, gamma_convert=False)

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
