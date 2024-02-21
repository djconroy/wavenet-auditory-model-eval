"""
This script is a modified version of models/wavenet/run_model.py,
which was written by Anil Nagathil.

Significant changes to the code are indicated below (insignificant changes aren't).
"""
import torch
import numpy as np
import yaml
import librosa
from scipy.io import savemat
from argparse import ArgumentParser
from models.wavenet.classes import WaveNet
from models.wavenet.utils import utils

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("--file",
                        type=str,
                        default='./models/wavenet/audio/61-70968-0000.flac',
                        help='Audio file to be processed')
    parser.add_argument("--spl",
                        type=float,
                        default=60,
                        help='Sound pressure level at which the signal is processed')
    parser.add_argument("--segdur",
                        type=float,
                        default=1,
                        help='Segment duration')
    parser.add_argument("--model",
                        type=str,
                        default="./models/wavenet/model/musan31rfa3-1fullSet_20231014-145738.pt",
                        help='Path to model parameters (.pt file)')
    parser.add_argument("--config",
                        type=str,
                        default="./models/wavenet/config/config31rfa3-1fullSet.yaml",
                        help='Path to config file (.yaml file)')
    parser.add_argument("--proc",
                        type=str,
                        choices=["cpu", "cuda"],
                        default="cpu",
                        help='Choose cpu or gpu processing')
    args = parser.parse_args()

    # load configuration file
    with open(args.config, 'r') as ymlfile:
        conf = yaml.safe_load(ymlfile)

    # constants
    sigMax = torch.tensor(55)
    ihcogramMax = torch.tensor(1.33)
    ihcogramMax = utils.comp(ihcogramMax, conf['scaleWeight'], conf['scaleType'])
    fs = 16000

    # number of samples to be skipped due to WaveNet processing
    skipLength = (2 ** conf['nLayers']) * conf['nStacks']

    # select processing device (either cpu or cuda)
    device = torch.device(args.proc)

    ## initialize WaveNet and load model paramaters
    NET = WaveNet.WaveNet(conf['nLayers'],
                          conf['nStacks'],
                          conf['nChannels'],
                          conf['nResChannels'],
                          conf['nSkipChannels'],
                          conf['numOutputLayers'])
    NET.to(device)
    NET.load_state_dict(torch.load(args.model, map_location=torch.device('cpu')))

    #####################################################
    ############ START OF CODE WRITTEN BY ME ############
    #####################################################

    print("Running run_wavenet_model.py")

    """
    Hack:
    With sr=None, librosa.load preserves the native sampling rate of the input,
    which is 16000 Hz for TIMIT audio signals, which is the expected sampling rate
    for inputs to the WaveNet model.
    """
    signal, _ = librosa.load(args.file, sr=None)

    """
    Hack:
    The speech signal in "noisy.wav" has already been normalized to an SPL,
    so only other speech signals need to be normalized to the specified SPL.
    """
    if args.file != "noisy.wav":
        rms = np.sqrt(np.mean(np.square(signal)))
        signal = (signal / rms) * 20e-6 * (10 ** (args.spl / 20))

    # Convert the signal from a numpy array to a tensor
    signal = torch.from_numpy(signal)

    #####################################################
    ############ END OF CODE WRITTEN BY ME ##############
    #####################################################

    signal = signal.to(device)

    # segmentation
    frame_shift = int(args.segdur * fs)
    frame_len = frame_shift + skipLength
    sigLen = signal.shape[0]

    ihcogram_pred = np.zeros((80, sigLen - skipLength + 1))

    start_idx = 0
    end_idx = 0
    frame_idx = 0

    # loop over frames
    while end_idx < sigLen:
        start_idx = frame_idx * frame_shift
        end_idx = np.min([sigLen, frame_idx * frame_shift + frame_len])

        if sigLen - end_idx < skipLength:
            end_idx = sigLen
        signal_seg = signal[start_idx:end_idx]
        signal_seg = signal_seg.to(device)

        # process signal segment and scale back to original range
        ihcogram_pred_seg = NET(signal_seg)
        ihcogram_pred_seg = ihcogram_pred_seg * ihcogramMax
        ihcogram_pred_seg = utils.invcomp(ihcogram_pred_seg, conf['scaleWeight'], conf['scaleType'])

        # write IHCogram segment into full IHCogram array
        ihcogram_pred_seg = ihcogram_pred_seg.cpu().detach().numpy()
        ihcogram_pred[:, start_idx:end_idx - skipLength + 1] = ihcogram_pred_seg

        frame_idx += 1

    # save IHCogram for MATLAB comparison
    savemat('tmp.mat', {'ihcogram': ihcogram_pred})
