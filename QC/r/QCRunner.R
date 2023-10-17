library(plumber)

path <- "/app/QCProcessing.R"
print(path)
r <- plumb(path)
r$run(port = 8001)