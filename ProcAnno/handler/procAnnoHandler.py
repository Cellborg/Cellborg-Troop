import requests
import boto3
import json
import os

client = boto3.client('sqs', region_name='us-west-2')
sns = boto3.client('sns', region_name='us-west-2')
queue_url = os.environ.get("SQS_QUEUE_URL")
env = os.environ.get("ENVIRONMENT", "dev")

print(f"Cellborg Troop - Processing and Annotation Python container running in {env} environment")

if env == "dev":
    SNS_TOPIC = 'arn:aws:sns:us-west-2:865984939637:QCCompleteTopic'
else:
    SNS_TOPIC = f'arn:aws:sns:us-west-2:865984939637:QCComplete-{env}-Topic'

def send_sns(data):
    topic_arn = SNS_TOPIC
    print(f'Sending {data} SNS message to {SNS_TOPIC}')
    response = sns.publish(
        TopicArn = topic_arn,
        Message = json.dumps(data),
    )
    return response

def send_request(endpoint, data):
    url = f"http://127.0.0.1:8002{endpoint}" 
    headers = {'Content-Type': 'application/json'}
    response = requests.post(url, json=data, headers=headers)
    return response.json()

def send_shutdown_request():
    url = f"http://127.0.0.1:8002/shutdown"
    try:
        response = requests.post(url, json = {"status": "complete"}, headers = {'Content-Type': 'application/json'})
        response.raise_for_status()
    except requests.ConnectionError:
        print("Connection was closed by the server (expected behavior during shutdown).")
    except requests.RequestException as e:
        print(f"An error occurred: {e}")

projectInitialized = False

while True:
    response = client.receive_message(QueueUrl=queue_url, MaxNumberOfMessages=10, WaitTimeSeconds=10, VisibilityTimeout=900)
    print(response)
    if 'Messages' not in response:
        print("No Message")
        continue

    for message in response['Messages']:

        try:
            queen_service_request_raw_data = message['Body']

            queen_service_request = json.loads(queen_service_request_raw_data)
            print(queen_service_request)
            request_type = queen_service_request['requestType']
            project = queen_service_request["project"]
            user = queen_service_request["user"]

            if request_type == "initializeProject":
                print("Initializing annData object")
                datasets = queen_service_request["datasets"]

                response = send_request('/init_endpoint', {"user": user, "project": project, "datasets": datasets})
                if response["success"]:
                    print("Initializing Project Successful... Sending SNS message")
                    project_initialized = True
                    data = {
                        "user": user, 
                        "project": project, 
                        "completed_step": "Initialize"
                    }
                    response = send_sns(data)
                    print(response)
            if request_type == "clustering":
                print("Beginning clustering now...")
                resolution = queen_service_request['resolution']

                response = send_request('/clustering', {"user":user, "project":project, "resolution":resolution})
                if response['success']:
                    print("Clustering was successful...")
                    data = {
                        "user":user,
                        "project":project,
                        "geneNames": response["gene_names"],
                        "clusters": response["clusters"],
                        "completed_step":"Clustered"
                    }
                    response = send_sns(data)
                    print(response)
                
            if request_type == "annotations":
                print("Beginning annotations now...")
                
            else:
                print("Invalid Request Type: ", request_type)
            


        except Exception as e:
            print("Error:", e)
        finally:
            client.delete_message(QueueUrl=queue_url, ReceiptHandle=message['ReceiptHandle'])


