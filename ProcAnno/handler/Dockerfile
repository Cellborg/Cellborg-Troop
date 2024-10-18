FROM python:3.12.1


#create /app dir
RUN mkdir /procHandler
#make working directory /app
WORKDIR /procHandler
COPY ./requirements.txt /procHandler/requirements.txt
RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

COPY ./procAnnoHandler.py /procHandler
CMD ["python", "/app/QCHandler.py"]

