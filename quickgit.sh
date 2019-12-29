echo "$1"
echo "Is the above commit message correct?"
options=("Type 1 for yes and enter." "Type 2 for no and enter.")
select yn in "${options[@]}"; do
  case $yn in
    "Type 1 for yes and enter.") 
      git init; git add .; git commit -m "$1"; git push origin master; break;;
    "Type 2 for no and enter.") break;;
  esac
done
