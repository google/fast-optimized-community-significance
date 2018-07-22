# Copyright 2018 Google LLC
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     https://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
 
context("FOCS Backend Check")
library(igraph)
RUN_SEED <- 123454321

test_that("A random community is scored appropriately", {
  set.seed(RUN_SEED)
  network <- degree.sequence.game(rep(5, 20), method="simple")
  membership <- cluster_louvain(network)$membership
  comm <- which(membership == 1)
  score <- unlist(FOCS(list(comm), network, p=0.25))
  expect_equal(signif(score, 7), 0.2162592)
})
