import tensorflow as tf

bitSize = 16

# 2 input vectors
x = tf.placeholder(tf.uint8, [None, 2])

weights = tf.Variable(tf.zeros[2, bitSize])
bias = tf.Variable(tf.zeros([bitSize]))

y = tf.nn.softmax(tf.matmul(x, weights) + b)

y_ = tf.placeholder(tf.uint8, [None, bitSize])
cross_entropy = tf.reduce_mean(-tf.reduce_sum(y_ * tf.log(y), reduction_indices=[1]))

train_step = tf.train.GraidientDescentOptimizer(0.5).minimize(cross_entropy)

init = tf.initialize_all_variables()

sess = tf.Session()
sess.run(init)

for i in range(1000):
  # batch_xs, batch_yz = mnist.train.next_batch(100)
  sess.run(train_step, feed_dict={x: batch_xs, y_: batch_yz})
