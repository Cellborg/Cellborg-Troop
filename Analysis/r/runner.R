library(plumber)

path <- "/app/queen.R"
print(path)
r <- plumb(path)
r$run(port = 8000)