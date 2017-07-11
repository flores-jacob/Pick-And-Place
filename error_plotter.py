import json
import numpy as np
from plotly.graph_objs import Box
from plotly.offline import plot

# from IK import T_total

# compute the FK, then compute the errors
# FK_matrix = T_total.evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4: theta4, q5: theta5, q6: theta6})

# how to jsonify numpy arrays https://stackoverflow.com/a/32850511
theta1_list = np.random.uniform(low = -0.1, high=0.1, size=50).tolist()
theta2_list = np.random.uniform(low = -0.1, high=0.1, size=50).tolist()
theta3_list = np.random.uniform(low = -0.1, high=0.1, size=50).tolist()
theta4_list = np.random.uniform(low = -0.1, high=0.1, size=50).tolist()
theta5_list = np.random.uniform(low = -0.1, high=0.1, size=50).tolist()
theta6_list = np.random.uniform(low = -0.1, high=0.1, size=50).tolist()

error_data = {
    "theta1": theta1_list,
    "theta2": theta2_list,
    "theta3": theta3_list,
    "theta4": theta4_list,
    "theta5": theta5_list,
    "theta6": theta6_list,
}

# save to disk
error_list_path = "./error_list.json"
with open(error_list_path, "w") as error_list_file:
    json_error_list = json.dump(error_data, error_list_file)

# open error list from disk
with open(error_list_path, "r") as error_list_file:
    json_error_list = json.loads(error_list_file.read())

print ("json error list", json_error_list)


theta1_list = json_error_list["theta1"]
theta2_list = json_error_list["theta2"]
theta3_list = json_error_list["theta3"]
theta4_list = json_error_list["theta4"]
theta5_list = json_error_list["theta5"]
theta6_list = json_error_list["theta6"]


# prepare the data for the box plot for plotly to plot
# box plot taken here https://plot.ly/python/box-plots/
theta1_errors = Box(x=1,
    y=theta1_list,
    name='theta1'
)
theta2_errors = Box(
    x=2,
    y=theta2_list,
    name='theta2'

)
theta3_errors = Box(
    x=3,
    y=theta3_list,
    name='theta3'

)
theta4_errors = Box(
    x=4,
    y=theta4_list,
    name='theta4'

)
theta5_errors = Box(
    x=5,
    y=theta5_list,
    name='theta5'

)
theta6_errors = Box(
    x=6,
    y=theta6_list,
    name='theta6'
)
data = [theta1_errors, theta2_errors, theta3_errors, theta4_errors, theta5_errors, theta6_errors]
plot({"data": data})