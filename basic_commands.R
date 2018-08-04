### Basic commands

library(workflowr)

# Find all the markdown files in the project which have been updated, and re-build them.
wflow_build()
# Find the status of the markdown files and whether they're ready to be published, or whether they are up-to-date
wflow_status()
# 'Unp = Unpublished' <- html files not committed to git repo yet
wflow_publish("analysis/*Rmd"
              ,message = "* Insert a commit message *"
              ,republish = T) # to publish the unpublished stuff

# We've only committed(published) to our local repository.
# Create a remote (and call it the default - origin)
# # Only need to run it once (at the start)
# wflow_git_remote(remote = "origin"
#                  ,user = "blaw36"
#                  ,repo = "Masters_Project")
# Now we'll push to our git repo online.
wflow_git_push()

# After setting up my SSH (to avoid continuous prompting of usernames and passwords):
wflow_git_remote(remote = "origin"
                 , user = "blaw36"
                 , repo = "Masters_Project",
                 protocol = "ssh", action = "set_url")
