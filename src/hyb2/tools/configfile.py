import argparse, yaml

def parse_with_config(parser, argv):

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=None)

    pre_args, _ = pre.parse_known_args(argv)
    if pre_args.config:
        with open(pre_args.config) as fh:
            cfg = yaml.safe_load(fh) or {}
        valid = {a.dest for a in parser._actions if a.dest not in ("help", "config")}
        unknown = set(cfg) - valid
        if unknown:
            parser.error(f"unknown config keys: {sorted(unknown)}")
        parser.set_defaults(**cfg)

        for action in parser._actions:
            if action.dest in cfg:
                action.required = False

    return parser.parse_args(argv)