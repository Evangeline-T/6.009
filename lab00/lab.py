# No Imports Allowed!


from git import safe_decode
from matplotlib.pyplot import sca


def backwards(sound):
    back = sound.copy()
    back['samples'] = back['samples'][::-1]
    return back


def mix(sound1, sound2, p):
    if sound1['rate'] != sound2['rate']:
        return None
    mixed = {'rate': sound1['rate'], 'samples': []}
    for s1, s2 in zip(sound1['samples'], sound2['samples']):
        mixed['samples'].append(p * s1 + (1 - p) * s2)
    return mixed


def convolve(sound, kernel):
    sound_len = len(sound['samples'])
    kernel_len = len(kernel)
    convolved = [ 0 for _ in range(sound_len + kernel_len - 1) ]
    for i in range(kernel_len):
        for j in range(sound_len):
            convolved[i + j] += kernel[i] * sound['samples'][j]
    return {'rate': sound['rate'], 'samples': convolved}


def echo(sound, num_echoes, delay, scale):
    sample_delay = round(delay * sound['rate'])
    sound_len = len(sound['samples'])
    sample_len =  num_echoes * sample_delay + sound_len
    new_sample = [ 0 for _ in range(sample_len) ]
    echos = sound['samples'].copy()
    for i in range(0, sample_len, sample_delay):
        if i + sound_len > sample_len:
            break
        for j in range(sound_len):
            new_sample[i + j] += echos[j]
        echos = [echo * scale for echo in echos]
    return {'rate': sound['rate'], 'samples': new_sample}

def echo_conv(sound, num_echoes, delay, scale):
    # warning: very slow
    sample_delay = round(delay * sound['rate'])
    echo_filter = [1] + sum([[*[0 for _ in range(sample_delay - 1)], scale ** i] for i in range(1, num_echoes + 1)], [])
    return convolve(sound, echo_filter)


def pan(sound):
    left = sound['left'][:]
    right = sound['right'][:]
    dur = len(left)
    for i in range(dur):
        left[i] *= (1 - i / (dur - 1))
        right[i] *= ( i / (dur - 1))
    return {'rate': sound['rate'], 'left': left, 'right': right}


def remove_vocals(sound):
    sample = [l - r for l, r in zip(sound['left'], sound['right'])]
    return {'rate': sound['rate'], 'samples': sample}


def bass_boost_kernel(N, scale=0):
    """
    Construct a kernel that acts as a bass-boost filter.

    We start by making a low-pass filter, whose frequency response is given by
    (1/2 + 1/2cos(Omega)) ^ N

    Then we scale that piece up and add a copy of the original signal back in.
    """
    # make this a fake "sound" so that we can use the convolve function
    base = {'rate': 0, 'samples': [0.25, 0.5, 0.25]}
    kernel = {'rate': 0, 'samples': [0.25, 0.5, 0.25]}
    for i in range(N):
        kernel = convolve(kernel, base['samples'])
    kernel = kernel['samples']

    # at this point, the kernel will be acting as a low-pass filter, so we
    # scale up the values by the given scale, and add in a value in the middle
    # to get a (delayed) copy of the original
    kernel = [i * scale for i in kernel]
    kernel[len(kernel)//2] += 1

    return kernel


# below are helper functions for converting back-and-forth between WAV files
# and our internal dictionary representation for sounds

import io
import wave
import struct

def load_wav(filename, stereo=False):
    """
    Given the filename of a WAV file, load the data from that file and return a
    Python dictionary representing that sound
    """
    f = wave.open(filename, 'r')
    chan, bd, sr, count, _, _ = f.getparams()

    assert bd == 2, "only 16-bit WAV files are supported"

    out = {'rate': sr}

    if stereo:
        left = []
        right = []
        for i in range(count):
            frame = f.readframes(1)
            if chan == 2:
                left.append(struct.unpack('<h', frame[:2])[0])
                right.append(struct.unpack('<h', frame[2:])[0])
            else:
                datum = struct.unpack('<h', frame)[0]
                left.append(datum)
                right.append(datum)

        out['left'] = [i/(2**15) for i in left]
        out['right'] = [i/(2**15) for i in right]
    else:
        samples = []
        for i in range(count):
            frame = f.readframes(1)
            if chan == 2:
                left = struct.unpack('<h', frame[:2])[0]
                right = struct.unpack('<h', frame[2:])[0]
                samples.append((left + right)/2)
            else:
                datum = struct.unpack('<h', frame)[0]
                samples.append(datum)

        out['samples'] = [i/(2**15) for i in samples]

    return out


def write_wav(sound, filename):
    """
    Given a dictionary representing a sound, and a filename, convert the given
    sound into WAV format and save it as a file with the given filename (which
    can then be opened by most audio players)
    """
    outfile = wave.open(filename, 'w')

    if 'samples' in sound:
        # mono file
        outfile.setparams((1, 2, sound['rate'], 0, 'NONE', 'not compressed'))
        out = [int(max(-1, min(1, v)) * (2**15-1)) for v in sound['samples']]
    else:
        # stereo
        outfile.setparams((2, 2, sound['rate'], 0, 'NONE', 'not compressed'))
        out = []
        for l, r in zip(sound['left'], sound['right']):
            l = int(max(-1, min(1, l)) * (2**15-1))
            r = int(max(-1, min(1, r)) * (2**15-1))
            out.append(l)
            out.append(r)

    outfile.writeframes(b''.join(struct.pack('<h', frame) for frame in out))
    outfile.close()


if __name__ == '__main__':
    # code in this block will only be run when you explicitly run your script,
    # and not when the tests are being run.  this is a good place to put your
    # code for generating and saving sounds, or any other code you write for
    # testing, etc.

    # here is an example of loading a file (note that this is specified as
    # sounds/hello.wav, rather than just as hello.wav, to account for the
    # sound files being in a different directory than this file)
    # mystery = load_wav('sounds/mystery.wav')
    # write_wav(backwards(mystery), 'mystery_reversed.wav')

    # synth = load_wav('sounds/synth.wav')
    # water = load_wav('sounds/water.wav')
    # write_wav(mix(synth, water, 0.2), 'mix_of_synth_water.wav')

    # ice_and_chilli = load_wav('sounds/ice_and_chilli.wav')
    # write_wav(convolve(ice_and_chilli, bass_boost_kernel(1000, 1.5)), 'convolved_ice.wav')
    
    # chord = load_wav('sounds/chord.wav')
    # write_wav(echo(chord, 5, 0.3, 0.6), 'echoed_chord.wav')

    # car = load_wav('sounds/car.wav', stereo=True)
    # write_wav(pan(car), 'car_pan.wav')

    # lookout_mountain = load_wav('sounds/lookout_mountain.wav', stereo=True)
    # write_wav(remove_vocals(lookout_mountain), 'lookout_mountain_remove_vocal.wav')
    pass