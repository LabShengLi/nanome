import logging

from global_config import data_base_dir

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)

    logging.info("Hello")
    logging.info(data_base_dir)
