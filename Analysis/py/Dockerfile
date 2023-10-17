FROM python:3.11-slim

RUN mkdir /app
WORKDIR /app

RUN pip install requests boto3

COPY knight.py /app/knight.py

CMD ["python", "/app/knight.py"]
