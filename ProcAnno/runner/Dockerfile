FROM python:3.12.1

WORKDIR /runner

COPY ./requirements.txt /runner/requirements.txt
RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

COPY ./app /runner/app

CMD ["fastapi", "run", "runner/main.py", "--port", "8001"]