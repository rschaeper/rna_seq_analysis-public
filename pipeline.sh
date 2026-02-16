#!/bin/bash

# define the name of the global functions scripts
GLOBAL_FUNCTIONS_SCRIPT="functions.R"

# check if the global functions script exists
if [ ! -f "$GLOBAL_FUNCTIONS_SCRIPT" ]; then
  echo "Error: Global functions script '$GLOBAL_FUNCTIONS_SCRIPT' not found in the current directory."
  echo "Please ensure it's in the same directory or provide the full path."
  exit 1
fi

# loop through R scripts in numerical order,
# explicitly excluding the global functions script itself and "unused_code.R",
# and ensuring only scripts starting with a number are processed.
for script in $(ls [0-9]*.R | sort | grep -v "$GLOBAL_FUNCTIONS_SCRIPT" | grep -v "unused_code.R"); do
  echo "Running $script..."
  # execute Rscript:
  # -e "source('$GLOBAL_FUNCTIONS_SCRIPT');" sources the global functions.
  # -e "source('$script')" then sources the current pipeline script.
  # both are executed in the same R session for this iteration.
  Rscript -e "source('$GLOBAL_FUNCTIONS_SCRIPT'); source('$script')"

  # check the exit status of the Rscript command
  if [ $? -ne 0 ]; then
    echo "Error in $script. Stopping execution."
    exit 1
  fi
done

echo "All scripts executed successfully."