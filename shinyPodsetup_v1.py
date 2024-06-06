import argparse
import glob
import os


def main():
    parser = argparse.ArgumentParser(
        description="Custom script with command-line parameters"
    )
    parser.add_argument(
        "-c", "--copy", required=True, help="Specify the LATCH URI to copy."
    )
    parser.add_argument(
        "-a",
        "--activate",
        action="store_true",
        help="Activate the command (without this, it's just a practice run)"
    )
    parser.add_argument(
        "-l", "--location", default="/root/", help="Optional location argument"
    )

    args = parser.parse_args()

    customAppPath = "/opt/latch/custom_app"
    if not args.activate:

        print(
            "\n\n**********\n*\tRunning in Dress Rehersal Mode.\n*\tFiles \
            WILL be copied.\n*\tApp will NOT be activated.\n*\tAdd the -a or \
            --activate to cmd line.\n**********\n"
            )

    # destination folder
    topCopiedDirName = args.copy.split("/")[-1]
    topCopiedPath = f"{args.location}/{topCopiedDirName}/".replace("//", "/")
    fileCopiedAlready = os.path.exists(topCopiedPath)
    if not fileCopiedAlready:

        # Execute the main command specified by -c
        cpCmd = f"latch cp {args.copy} {args.location}"
        print(f"Copying {args.copy} to destination {topCopiedPath}\n")
        print(f"Running main command:\n\t{cpCmd}")
        os.system(cpCmd)
        print(">>DONE\n")

    else:

        print(
            f"\n\n!!Directory already exists [{topCopiedPath}]. Remove if you \
            want to re-copy.\n"
        )

    print(f"Searching in {topCopiedPath} for 'ui.R' file")
    files = glob.glob(f"{topCopiedPath}**/ui.R", recursive=True)

    if len(files) == 1:

        shinyDir = os.path.dirname(files[0]) + "/"
        print(f"\tShinyApp Directory Found: {shinyDir}\n")

        serviceFile = (
            "#!/usr/bin/env Rscript\n\n"
            "#Need to specify lib that Rstudio is using as it's different "
            "from R on commandline\n"
            ".libPaths('/usr/local/lib/R/site-library');"
            f"shiny::runApp('{shinyDir}', "
            "port=5000, host='0.0.0.0')\n"
        )

        tmpFilePath = "/root/tmp_custom_app.txt"
        print(f"Writing new file [{tmpFilePath}] to /root/tmp_custom_app.txt")

        with open(tmpFilePath, mode="w") as file:
            file.write(serviceFile)

        print(">>DONE\n")
        print(
            f"Backing up existing file [{customAppPath}] to \
            /root/custom_app.old.txt"
        )
        os.system(f"cp {customAppPath} /root/custom_app.old.txt")
        print(">>DONE\n")

        if args.activate:

            print(f"Copying new file to {customAppPath}")
            os.system(f"cp /root/tmp_custom_app.txt {customAppPath}")
            print(">>DONE\n")

            actCmd = "systemctl --user enable latch-custom-app"
            startCmd = "systemctl --user start latch-custom-app"

            print(f"Enabling App to start upon pod startup: {actCmd}")
            os.system(actCmd)
            print(">>DONE\n")

            print(f"Starting up cusom_app service: {startCmd}")
            os.system(startCmd)
            print(">>DONE\n")

            print(
                "!! You can check to see if service is active by running: \
                >systemctl --user status latch-custom-app on the command line"
            )
        else:
            print(
                "****App was not activated. Run with -a or --activate on \
                commandline.****\n\n"
            )
    else:
        print(
            "Could not find obvious ShinyApp directory containing the 'ui.R' \
            file"
        )
        for file in files:
            print(f"Found: {file}")


if __name__ == "__main__":
    main()
