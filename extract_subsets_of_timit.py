import os
import shutil

TIMIT_DATASETS_DIR = "datasets/timit/"
TIMIT_CORE_TEST_SET = "core-test-set/"
TIMIT_CORE_TEST_SET_PREFIX = "TIMIT/TIMIT/TEST/DR"
TIMIT_CORE_TEST_SET_SPEAKERS = {
    "F": ["ELC0", "PAS0", "PKT0", "JLM0", "NLP0", "MGD0", "DHC0", "MLD0"],
    "M1": ["DAB0", "TAS1", "JMP0", "LLL0", "BPM0", "CMJ0", "GRT0", "JLN0"],
    "M2": ["WBT0", "WEW0", "LNT0", "TLS0", "KLT0", "JDH0", "NJM0", "PAM0"],
}
TIMIT_CORE_TEST_SET_DIRS = {
    key: [
        f"{TIMIT_CORE_TEST_SET_PREFIX}{dialect_region}/{key[0]}{speaker}"
        for dialect_region, speaker in enumerate(speakers, start=1)
    ]
    for key, speakers in TIMIT_CORE_TEST_SET_SPEAKERS.items()
}

for key, core_test_set in TIMIT_CORE_TEST_SET_DIRS.items():

    core_test_set_dir = TIMIT_DATASETS_DIR + TIMIT_CORE_TEST_SET + key

    for dialect_region, speaker in enumerate(TIMIT_CORE_TEST_SET_SPEAKERS[key], start=1):

        # Create a directory to store the test data files for the speaker
        core_test_set_speaker_dir = f"{core_test_set_dir}/DR{dialect_region}/{speaker}"
        os.makedirs(core_test_set_speaker_dir, exist_ok=True)

        # Get the test data files for the speaker
        timit_data = os.listdir(core_test_set[dialect_region - 1])
        # Exclude the SA sentence data files
        timit_data = [file for file in timit_data if not file.startswith("SA")]

        # Copy the test data files for the speaker into the newly-created directory
        for file in timit_data:
            #shutil.move(f"{core_test_set[dialect_region - 1]}/{file}", core_test_set_speaker_dir)
            shutil.copy(f"{core_test_set[dialect_region - 1]}/{file}", core_test_set_speaker_dir)
