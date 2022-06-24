import matplotlib.pyplot as plt
import pickle

with open('para_changing_z_feature', 'rb') as filename:
  s_x_test_z = pickle.load(filename)

with open('output_field_az_changing_z_feature', 'rb') as filename:
  s_y_test_z = pickle.load(filename)

Z = []
for elem in s_x_test_z:
  Z.append(elem[5])

print ('With all other features keeping constant: ')
print ('Changing Z feature with a fix step.')
for i in range(5):
    print('With Z value of: ' + str(Z[i]) + ' has M.field of: ' + str(s_y_test_z[i]))

plt.scatter(Z, s_y_test_z, color = 'blue')
plt.title('Magnetic field with changing Z-feature')
plt.xlabel('Z')
plt.show()