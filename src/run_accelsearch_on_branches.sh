#!/bin/bash

original_dir=$(pwd)

presto_dir="presto/src"
branches=("master" "cufftPlanManyCache" "2132b2100427963564050b1b27dc1ff4c5c233b1" "feature/20240219/dynamicBatch")
cmd_params=(
  "-zmax 20 -wmax 20"
  "-zmax 50 -wmax 50"
  "-zmax 100 -wmax 100"
  "-zmax 200 -wmax 200"
)


for branch in "${branches[@]}"
do
  cd $presto_dir || exit
  echo "Checking out branch $branch..."
  git checkout $branch > /dev/null 2>&1
  echo "Building on branch $branch, output is suppressed..."
  make cleaner > /dev/null 2>&1 && make accelsearch_cu accelsearch > /dev/null 2>&1
  cd "$original_dir" || exit
  echo "Finished building on branch $branch."

  for params in "${cmd_params[@]}"
  do
    echo "Running with params: $params"
    total_time=0
    valid_runs=0 # Track number of successful runs
    for i in {1..3}
    do
      output=$(accelsearch_cu $params FRB121102_tracking-M01_0706_ds1_0_DM99.70.fft 2>&1)
      time=$(echo "$output" | grep -oP 'Total time: \K[0-9.]+')
      if [ -z "$time" ]; then
        echo "No time found for run $i with params: $params"
      else
        total_time=$(echo "$total_time + $time" | bc)
        valid_runs=$((valid_runs + 1))
      fi
    done

    if [ $valid_runs -eq 0 ]; then
      echo "No valid runs completed for params $params."
    else
      average_time=$(echo "scale=3; $total_time / $valid_runs" | bc)
      echo "Average Time for params $params: $average_time sec"
    fi
  done
done
