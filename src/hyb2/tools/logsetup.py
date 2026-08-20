import logging

def configure_logging(verbose: bool = False):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(message)s"
    )

    for noisy in ("numba",):
        logging.getLogger(noisy).setLevel(logging.WARNING)

def is_quiet() -> bool:
    # True when DEBUG logging is off (i.e. not -V)
    return not logging.getLogger().isEnabledFor(logging.DEBUG)