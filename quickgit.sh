echo "$1"
echo "Is the above commit message correct?"
options=("y" "n")
select yn in "${options[@]}"; do
  case $yn in
    y) git init; git add .; git commit -m "$1"; git push origin master; break;;
    n) break;;
  esac
done
