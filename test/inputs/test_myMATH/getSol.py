#%%
import numpy as np

# %%
x_list = []
y_list = []
coeff_list = []
order_list = []

# linear function
x = np.linspace(0, 1, 2)
y = 4*x + 1
x_list.append(x)
y_list.append(y)
coeff_list.append(np.polyfit(x, y, 1)[::-1])
order_list.append(1)

x = np.linspace(0, 1, 10)
y = 2.5*x + 10
x_list.append(x)
y_list.append(y)
coeff_list.append(np.polyfit(x, y, 1)[::-1])
order_list.append(1)

# quadratic function
x = np.linspace(0, 1, 10)
y = 5.5*x**2 + 4*x + 1
x_list.append(x)
y_list.append(y)
coeff_list.append(np.polyfit(x, y, 2)[::-1])
order_list.append(2)

# cubic function
x = np.linspace(0, 1, 10)
y = 3.5*x**3 + 5.5*x**2 + 4*x + 1
x_list.append(x)
y_list.append(y)
coeff_list.append(np.polyfit(x, y, 3)[::-1])
order_list.append(3)

# with noise 
x = np.linspace(0, 1, 10)
y = 3.5*x**3 + 5.5*x**2 + 4*x + 1 + np.random.normal(0, 1e-3, 10)
x_list.append(x)
y_list.append(y)
coeff_list.append(np.polyfit(x, y, 3)[::-1])
order_list.append(3)


#%%
with open("x_list.txt", "w") as f:
    for x in x_list:
        f.write(",".join(map(str, x)) + '\n')

with open("y_list.txt", "w") as f:
    for y in y_list:
        f.write(",".join(map(str, y)) + '\n')
        
with open("order_list.txt", "w") as f:
    for order in order_list:
        f.write(str(order) + '\n')

with open("../../solutions/test_myMATH/coeff_list.txt", "w") as f:
    for coeff in coeff_list:
        f.write(",".join(map(str, coeff)) + '\n')
        
# %%
