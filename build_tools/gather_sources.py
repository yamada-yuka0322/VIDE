from git import Repo

repo = Repo(".")

assert repo.bare == False

t = repo.tree()


for entry in t.travers():
  print entry
