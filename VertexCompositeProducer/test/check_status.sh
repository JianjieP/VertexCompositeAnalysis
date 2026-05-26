for d in crab_projects/crab_P*; do
  echo "Checking $d"
  crab status -d "$d"
done

