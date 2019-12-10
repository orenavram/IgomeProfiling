def verify_file_is_not_empty(file_path):
    import logging
    logger = logging.getLogger('main')
    # make sure that there are results and the file is not empty
    with open(file_path) as f:
        if len(f.read(10).strip()) == 0:
            # TODO: write error to a global error file
            msg = f'Input file is empty {file_path}'
            logger.error(msg)
            raise RuntimeError(msg)
