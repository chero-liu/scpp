import argparse

from scpp.tools import utils
from scpp.__init__ import __VERSION__, ASSAY_LIST


# Define a custom argument formatter class
class ArgFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


def main():
    """Main function for handling command-line arguments and subcommands"""
    parser = argparse.ArgumentParser(
        description="A Python library and command-line tool for single-cell multi-sample analysis",
        formatter_class=ArgFormatter,
    )

    # Add a version argument
    parser.add_argument("-v", "--version", action="version", version=__VERSION__)

    # Create subparsers for different assays
    subparsers = parser.add_subparsers(dest="subparser_assay")

    # Loop through available assays
    for assay in ASSAY_LIST:
        # Import the initialization module for the specific assay
        init_module = utils.find_assay_init(assay)

        # Create a subparser for the assay with a description from the init module
        subparser_1st = subparsers.add_parser(assay, description=init_module.__PREFIX__)

        # Add a second level of subparsers for assay steps
        subparser_2nd = subparser_1st.add_subparsers()

        for step in init_module.__STEPS__:
            # Import the function and options module for the assay step
            step_module = utils.find_step_module(assay, step)
            func = getattr(step_module, step)
            func_opts = getattr(step_module, f"get_opts_{step}")

            # Create a subparser for the assay step with a custom argument formatter
            parser_step = subparser_2nd.add_parser(step, formatter_class=ArgFormatter)

            # Add options for the assay step
            func_opts(parser_step, sub_program=True)

            # Set the default function to be called for this subcommand
            parser_step.set_defaults(func=func)
    # Parse command-line arguments
    args = parser.parse_args()

    # Check if no arguments or subcommands were given, and display help
    if len(args.__dict__) <= 1:
        # No arguments or subcommands were given.
        parser.print_help()
        parser.exit()
    else:
        # Call the appropriate function for the selected subcommand
        args.func(args)


if __name__ == "__main__":
    main()
