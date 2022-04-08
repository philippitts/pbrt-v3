"""Testing suite for time-of-flight tracer film luminosity output"""

import luminance


def main():
    # exr = luminance.load('C:/git/pbrt-v3/build/Debug/output/corner_directlighting.exr')
    # dat = luminance.load('C:/git/pbrt-v3/build/Debug/output/corner_groundtruth.dat')

    exr = luminance.load(('C:/git/pbrt-v3/build/Debug/output/corner_path.exr'))
    dat = luminance.load(('C:/git/pbrt-v3/build/Debug/output/corner_pathtof.dat'))

    luminance.compare(exr, dat)


if __name__ == "__main__":
    main()
