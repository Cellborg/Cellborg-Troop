import requests
import boto3
import json
import os

client = boto3.client('sqs', region_name='us-west-2')
sns = boto3.client('sns', region_name='us-west-2')
queue_url = os.environ.get("SQS_QUEUE_URL")
env = os.environ.get("ENVIRONMENT", "dev") 
print(f"Cellborg Troop - QC Python container running in {env} environment")
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
    url = f"http://127.0.0.1:8001{endpoint}"  # Replace with the actual IP address of container host
    response = requests.post(url, json = data, headers = {'Content-Type': 'application/json'})
    return response.json()

def send_shutdown_request():
    url = f"http://127.0.0.1:8001/shutdown"
    try:
        response = requests.post(url, json = {"status": "complete"}, headers = {'Content-Type': 'application/json'})
        response.raise_for_status()
    except requests.ConnectionError:
        print("Connection was closed by the server (expected behavior during shutdown).")
    except requests.RequestException as e:
        print(f"An error occurred: {e}")
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
            request_type = queen_service_request["requestType"]
            project = queen_service_request["project"]
            user = queen_service_request["user"]
            dataset = queen_service_request["dataset"]

            if request_type == "qualityControl":

                subset_min = queen_service_request["min"]
                subset_max = queen_service_request["max"]
                subset_mt = queen_service_request["mt"]
                qc_request = {
                    "user": user, 
                    "project": project, 
                    "dataset": dataset,
                    "min": subset_min,
                    "max": subset_max,
                    "mt": subset_mt
                }            
                print("Sending QC request...")
                response = send_request('/qc_endpoint', qc_request)
                if response.get("success"):
                    print("QC Successful... Sending SNS message to clear dataset as completed...")
                    cell_count = response.get("cell_count")
                    gene_count = response.get("gene_count")
                    # Send SNS message
                    data = {
                        "user": user, 
                        "project": project, 
                        "dataset": dataset,
                        "complete": True,
                        "cell_count": cell_count,
                        "gene_count": gene_count
                    } 
                    response = send_sns(data)
                    print(response)
                else:
                    print(f"Error in QC: {response.get('message')}")
            
            elif request_type == "killServer":

                print("Shutting down the R QC server...")
                send_shutdown_request()
                print("Shutting down the python handler")
                exit(0)

        except Exception as e:
            print("Error:", e)
        finally:
            client.delete_message(QueueUrl=queue_url, ReceiptHandle=message['ReceiptHandle'])