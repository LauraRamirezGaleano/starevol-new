#!/bin/bash

# SYNTAX:
# bash screate_classical_Laura.sh <mass> <z> <starevol_*> <model_name>
# z is metallicity (digits after the decimal separator, i.e., 012345 for a metallicity of 0.012345)
#
# EXAMPLES:
# bash screate_classical_Laura.sh 0.78 000123 classical TEST_LAURA
# bash screate_classical.sh 1.0 014227 classical solar_model
# bash screate_classical.sh 1.0 014227 solar_model

INSTDIR="/Users/laura/CALCULS"
DIRGRIDS="$INSTDIR/starevol/PHYSICS/NUCLEAR/NETGEN/EVOL"
DIRBATCH="$INSTDIR/starevol/BATCH"
DIRM="$INSTDIR/MODELS"
DIRV="$INSTDIR/RESEVOL"
DIRS="$INSTDIR/RESULTS"
#DIRW="$INSTDIR/starevol/PHYSICS/ATMOSPHERE/PHOENIX/" 
DIRMOD="$INSTDIR/starevol/PHYSICS/INIT"
DIRCREATE="$INSTDIR/starevol/SCREATE"
DIRCODE="$INSTDIR/starevol/CODE"

# -------------------
# GETTING BASIC DATA
# -------------------

# Ensure that all required arguments are provided
if [ -z "$2" ]; then
  echo "Enter mass (5.0 for a 5 sm star), metallicity (004 for Z=0.004), version (2.51...), and other (mloss...)"
  exit 1
fi

# Assign input arguments to variables
m=$1
Z=$2
v=$3
b=$4
ms=$m

# Check if an optional fifth argument is provided, and assign it to 'ms' if it exists
if [ -n "$5" ]; then
  ms=$5
fi

# Check for bad input in Z, specifically if it starts with "0."
if [[ "$Z" =~ ^0\. ]]; then
  echo "Bad input for Z = 0.$Z"
  exit 1
fi

# Format 'ms' to two decimal places and use it to create 'name' and 'dirname'
mstring=$(printf "%.2f" "$ms")
name="m${mstring}z${Z}_${v}_${b}"
dirname="M${mstring}Z${Z}_${v}_${b}"
echo "Directory: M${mstring}Z${Z}"

# --------------------
# CHECKING FOR ERRORS
# --------------------

# Adjusting the string format based on the mass
if (( $(echo "$m >= 10" | bc -l) )); then
  nb=2
  na=5
else
  nb=1
  na=4
fi

# Function to adjust string length, adding leading zeros if necessary
stringadjust() {
  local input=$1
  local left=$2
  local right=$3
  printf "%0${left}.${right}f" "$input"
}

# Adjust the format of 'm' and store it in 'mprint'
mprint=$(stringadjust "$m" "$nb" "$na")

# Check if required files and directories exist
if [ ! -e "$DIRCREATE/starevol_${v}.par" ]; then
  echo "$DIRCREATE/starevol_X.XX.par does not exist"
  exit 1
fi

if [ -d "$DIRM/$dirname" ]; then
  echo "$DIRM/$dirname already exists"
  exit 1
fi

# Additional formatted values for 'mtini', 'znew', and 'zznew'
mtini=$(stringadjust "$m" 2 6)
znew=$(stringadjust "0.$Z" 1 9)
zznew=$(stringadjust "0.$Z" 1 8)

