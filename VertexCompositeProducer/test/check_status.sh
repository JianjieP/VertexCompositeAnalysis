for d in crab_projects*/crab_*; do
  echo "Checking $d"
  crab status -d "$d"
done

