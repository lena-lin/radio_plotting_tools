import imageio
from glob import glob
import datetime
import click


def create_gif(filenames, duration, output_dir):
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
    output_file = '{}/Gif-{}.gif'.format(output_dir, datetime.datetime.now().strftime('%Y-%M-%d-%H-%M-%S'))
    imageio.mimsave(output_file, images, duration=duration)


@click.command()
@click.argument('input_card', type=click.Path(file_okay=True, dir_okay=False))
@click.option('-o', '--output_directory', type=click.Path(file_okay=False, dir_okay=True), default='.')
def main(
    input_card,
    output_directory,
        ):
    image_filenames = sorted(glob(input_card))
    # from IPython import embed; embed()
    #dates = image_filenames.split('_3C84_')[1].split('.')[0]
    #from IPython import embed; embed()
    create_gif(image_filenames, 1, output_directory)


if __name__ == '__main__':
    main()
