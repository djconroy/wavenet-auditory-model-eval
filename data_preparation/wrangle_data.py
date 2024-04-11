from pathlib import Path
import scipy.io
from itertools import chain

DATASET = "timit/core-test-set"

parent_dir = Path.cwd().parent

data_dir = Path(parent_dir, "outputs/datasets", DATASET)
comparison_metrics_files = data_dir.glob("*/*/*/*/*/*/comparison_metrics.mat")

comparison_metrics_dir = Path(parent_dir, "comparison_metrics", DATASET)
comparison_metrics_dir.mkdir(parents=True, exist_ok=True)

sentences_file = Path(comparison_metrics_dir, "sentences.csv")
words_file = Path(comparison_metrics_dir, "words.csv")
phonemes_file = Path(comparison_metrics_dir, "phonemes.csv")

with (sentences_file.open("w") as sentences_csv,
      words_file.open("w") as words_csv,
      phonemes_file.open("w") as phonemes_csv):

    first_columns = ["SPL", "SNR", "Sex", "Dialect Region", "Speaker ID", "Sentence ID"]
    sentence_columns = ["Sentence"]
    word_columns = ["Word", "Word Number"]
    phoneme_columns = ["Phoneme", "Phoneme Number"]
    last_columns = ["Mean TFS NSIM", "Mean ENV NSIM", "Mean SDR"]
    last_columns.extend("TFS NSIM " + str(num) for num in range(1, 77))
    last_columns.extend("ENV NSIM " + str(num) for num in range(1, 77))
    last_columns.extend("SDR " + str(num) for num in range(1, 81))

    sentences_csv.write(",".join(chain(first_columns, sentence_columns, last_columns)) + "\n")
    words_csv.write(",".join(chain(first_columns, word_columns, last_columns)) + "\n")
    phonemes_csv.write(",".join(chain(first_columns, phoneme_columns, last_columns)) + "\n")

    for comparison_metrics_file in comparison_metrics_files:

        # Example storage structure for comparison metrics files:
        # F/DR1/ELC0/SI756/SPL60/SNR-15/comparison_metrics.mat
        parts = comparison_metrics_file.parts
        snr = parts[-2][3:]
        spl = parts[-3][3:]
        sentence_ID = parts[-4]
        speaker_ID = parts[-5]
        dialect_region = parts[-6][2:]
        sex = parts[-7][0]
        first_entries = [spl, snr, sex, dialect_region, speaker_ID, sentence_ID]

        comparison_metrics = scipy.io.loadmat(comparison_metrics_file, squeeze_me=True)

        words = comparison_metrics["words"]
        sentence = " ".join(words)
        phonemes = comparison_metrics["phonemes"]

        SDRs = comparison_metrics["SDRs"]
        mean_SDR = comparison_metrics["mean_SDR"]
        env_NSIMs = comparison_metrics["env_NSIMs"]
        mean_env_NSIM = comparison_metrics["mean_env_NSIM"]
        tfs_NSIMs = comparison_metrics["tfs_NSIMs"]
        mean_tfs_NSIM = comparison_metrics["mean_tfs_NSIM"]

        sentences_csv.write(",".join(chain(
            first_entries,
            [sentence, str(mean_tfs_NSIM), str(mean_env_NSIM), str(mean_SDR)],
            [str(tfs_NSIM) for tfs_NSIM in tfs_NSIMs],
            [str(env_NSIM) for env_NSIM in env_NSIMs],
            [str(sdr) for sdr in SDRs]
        )) + "\n")

        words_SDRs = comparison_metrics["words_SDRs"].T
        words_mean_SDRs = comparison_metrics["words_mean_SDRs"]
        words_env_NSIMs = comparison_metrics["words_env_NSIMs"].T
        words_mean_env_NSIMs = comparison_metrics["words_mean_env_NSIMs"]
        words_tfs_NSIMs = comparison_metrics["words_tfs_NSIMs"].T
        words_mean_tfs_NSIMs = comparison_metrics["words_mean_tfs_NSIMs"]

        for word_num, word in enumerate(words):
            words_csv.write(",".join(chain(
                first_entries,
                [word, str(word_num + 1), str(words_mean_tfs_NSIMs[word_num]),
                 str(words_mean_env_NSIMs[word_num]), str(words_mean_SDRs[word_num])],
                [str(word_tfs_NSIM) for word_tfs_NSIM in words_tfs_NSIMs[word_num]],
                [str(word_env_NSIM) for word_env_NSIM in words_env_NSIMs[word_num]],
                [str(word_SDR) for word_SDR in words_SDRs[word_num]]
            )) + "\n")

        phonemes_SDRs = comparison_metrics["phonemes_SDRs"].T
        phonemes_mean_SDRs = comparison_metrics["phonemes_mean_SDRs"]
        phonemes_env_NSIMs = comparison_metrics["phonemes_env_NSIMs"].T
        phonemes_mean_env_NSIMs = comparison_metrics["phonemes_mean_env_NSIMs"]
        phonemes_tfs_NSIMs = comparison_metrics["phonemes_tfs_NSIMs"].T
        phonemes_mean_tfs_NSIMs = comparison_metrics["phonemes_mean_tfs_NSIMs"]

        for phoneme_num, phoneme in enumerate(phonemes):
            phonemes_csv.write(",".join(chain(
                first_entries,
                [phoneme, str(phoneme_num + 1), str(phonemes_mean_tfs_NSIMs[phoneme_num]),
                 str(phonemes_mean_env_NSIMs[phoneme_num]), str(phonemes_mean_SDRs[phoneme_num])],
                [str(phoneme_tfs_NSIM) for phoneme_tfs_NSIM in phonemes_tfs_NSIMs[phoneme_num]],
                [str(phoneme_env_NSIM) for phoneme_env_NSIM in phonemes_env_NSIMs[phoneme_num]],
                [str(phoneme_SDR) for phoneme_SDR in phonemes_SDRs[phoneme_num]]
            )) + "\n")
