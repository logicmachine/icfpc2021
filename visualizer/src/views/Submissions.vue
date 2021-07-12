<template>
  <v-container>
    <v-row>
      <v-col>
        <h2 class="headline mb-3">Problem {{ problem_id }}</h2>
        <v-data-table :headers="headers" :items="submissions" class="elevation-1" :options="options">
          <template #item.id="{ item }">
            <router-link :to="{ name: 'Submission', params: { problem_id: problem_id, submission_id: item.id } }">
              {{ item.id }}
            </router-link>
          </template>
        </v-data-table>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
import axios from 'axios'

const api_options = {
  headers: {
    'Authorization': `Bearer ${process.env.VUE_APP_API_TOKEN}`
  }
}

export default {
  name: 'Submissions',

  data: () => {
    return {
      problem_id: null,
      headers: [
        { text: 'Submission', value: 'id' },
        { text: 'Dislikes', value: 'dislikes' },
        { text: 'Solver', value: 'solver' }
      ],
      submissions: [],
      options: {
        itemsPerPage: -1,
        sortBy: ['dislikes']
      }
    }
  },

  created(){
    this.fetchData()
  },

  watch: {
    $route: 'fetchData'
  },

  methods: {
    fetchData(){
      this.problem_id = this.$route.params.problem_id
      axios
        .get(`/api/submissions/${this.problem_id}`, api_options)
        .then(response => {
          this.submissions = response.data
        })
    }
  }
}
</script>
