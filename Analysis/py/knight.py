import requests
import boto3
import json
import os

client = boto3.client('sqs', region_name='us-west-2')
sns = boto3.client('sns', region_name='us-west-2')

queue_url = os.environ.get("SQS_QUEUE_URL")
env = os.environ.get("ENVIRONMENT", "dev") 
print(f"Cellborg Troop - Analysis Python container running in {env} environment")
if env == "dev":
    SNS_TOPIC = 'arn:aws:sns:us-west-2:865984939637:AnalysisStepCompleteTopic'
else:
    SNS_TOPIC = f'arn:aws:sns:us-west-2:865984939637:AnalysisStepComplete-{env}-Topic'

def send_sns(data):
    topic_arn = SNS_TOPIC
    print(f'Sending {data} SNS message to {SNS_TOPIC}')
    response = sns.publish(
        TopicArn = topic_arn,
        Message = json.dumps(data),
    )
    return response

def send_request(endpoint, data):
    url = f"http://127.0.0.1:8000{endpoint}" 
    headers = {'Content-Type': 'application/json'}
    response = requests.post(url, json=data, headers=headers)
    return response.json()

def send_shutdown_request():
    url = f"http://127.0.0.1:8000/shutdown"
    try:
        response = requests.post(url, json = {"status": "complete"}, headers = {'Content-Type': 'application/json'})
        response.raise_for_status()
    except requests.ConnectionError:
        print("Connection was closed by the server (expected behavior during shutdown).")
    except requests.RequestException as e:
        print(f"An error occurred: {e}")

project_initialized = False

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
            analysis_id = queen_service_request["analysisId"]

            if request_type == "initializeProject":
                print("Initializing seurat object")
                datasets = queen_service_request["datasets"]
                response = send_request('/init_endpoint', {"user": user, "project": project, "datasets": datasets, "analysis": analysis_id})
                if response.get("success"):
                    print("Initializing Project Successful... Sending SNS message")
                    project_initialized = True
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "Initialize"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "variableFeatures" and project_initialized:
                nFeatures = queen_service_request["nFeatures"]
                print("Finding variable features")
                response = send_request('/variableFeatures', {"user": user, "project": project, "analysis": analysis_id, "nFeatures": nFeatures})
                if response.get("success"):
                    print("Finding variable features was successful... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "VariableFeatures"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "cluster" and project_initialized:
                clustering_request_data = {
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id,
                    "cluster_neighbors": queen_service_request['neighbors'],
                    "cluster_clusters": queen_service_request['clusters'],
                    "cluster_dimensions": queen_service_request['dimensions'],
                    "cluster_reduction": queen_service_request['reduction']
                }
                print("Performing clustering")
                response = send_request('/cluster', clustering_request_data)
                if response.get("success"):
                    print("Performing clustering was successful... Sending SNS message") 
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "clusters"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "featurePlot" and project_initialized:
                gene_name = queen_service_request["gene_name"]
                feature_plot_request = {
                    "gene_name": gene_name,
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id
                }
                print("Finding gene feature plot", feature_plot_request)
                response = send_request('/genefeature', feature_plot_request)
                if response.get("success"):
                    print("Finding gene feature plot was successful... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "FeaturePlot",
                        "gene_name": gene_name
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "heatmap" and project_initialized:
                gene_names = queen_service_request["geneNames"]
                print(gene_names)
                heatmap_req = {
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id,
                    "geneNames": gene_names
                }
                response = send_request('/heatmap_analysis', heatmap_req)
                if response.get("success"):
                    print("Performing heatmap analysis was successful... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "Heatmap"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "psuedotime" and project_initialized:
                points = queen_service_request["points"]
                print(points)
                psuedotime_req = {
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id,
                    "points": points
                }
                response = send_request('/psuedotime', psuedotime_req)
                if response.get("success"):
                    print("Performing psuedotime analysis was successful... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "Psuedotime"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "dotplot" and project_initialized:
                genes = queen_service_request["genes"]
                print(genes)
                dotplot_req = {
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id,
                    "genes": genes
                }
                response = send_request('/dotplot', dotplot_req)
                if response.get("success"):
                    print("Collecting dot plot data was successful... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "DotPlot"
                    }
                    response = send_sns(data)
                    print(response)

            
            elif request_type == "vlnplots" and project_initialized:
                genes = queen_service_request["genes"]
                print(genes)
                vlnplots_req = {
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id,
                    "genes": genes
                }
                response = send_request('/violinplots', vlnplots_req)
                if response.get("success"):
                    print("Collecting violin plot data was successful... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "ViolinPlot"
                    }
                    response = send_sns(data)
                    print(response)
            
            elif request_type == "annotations" and project_initialized:
                annotations = queen_service_request["annotations"]
                print(annotations)
                annotations_req = {
                    "annotations": annotations
                }
                response = send_request('/annotations', annotations_req)
                if response.get("success"):
                    print("Annotating the data was successful... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "Annotations"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "allMarkers" and project_initialized:

                allMarkers_req = {
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id
                }
                response = send_request('/allMarkers', allMarkers_req)
                if response.get("success"):
                    print("Gathered all markers data... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "allMarkers"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "findMarkers" and project_initialized:
                cluster1 = queen_service_request["cluster1"]
                cluster2 = queen_service_request["cluster2"]
                print(cluster1)
                print(cluster2)
                findMarkers_req = {
                    "user": user, 
                    "project": project, 
                    "analysis": analysis_id,
                    "cluster1": cluster1,
                    "cluster2": cluster2
                }
                response = send_request('/findMarkers', findMarkers_req)
                if response.get("success"):
                    print("Gathered find markers data... Sending SNS message")
                    data = {
                        "user": user, 
                        "project": project, 
                        "analysisId": analysis_id,
                        "completed_step": "findMarkers"
                    }
                    response = send_sns(data)
                    print(response)

            elif request_type == "killServer":

                print("Shutting down the R QC server...")
                send_shutdown_request()
                print("Shutting down the python handler")
                exit(0)
                
            else:
                print("Invalid Request Type: ", request_type)
        except Exception as e:
            print("Error:", e)
        finally:
            client.delete_message(QueueUrl=queue_url, ReceiptHandle=message['ReceiptHandle'])