# -*- coding: UTF-8 -*-

import requests
from urllib.parse import urlencode
from pyquery import PyQuery as pq

host = 'm.weibo.cn'
base_url = 'https://%s/api/container/getIndex?' % host
user_agent = 'User-Agent: Mozilla/5.0 (iPhone; CPU iPhone OS 9_1 like Mac OS X) AppleWebKit/601.1.46 (KHTML, like Gecko) Version/9.0 Mobile/13B143 Safari/601.1 wechatdevtools/0.7.0 MicroMessenger/6.3.9 Language/zh_CN webview/0'

headers = {
    'Host': host,
    'Referer': 'https://m.weibo.cn/u/1665372775',
    'User-Agent': user_agent
}


# 按页数抓取数据
def get_single_page(page):
    params = {
        'type': 'uid',
        'value': 1665372775,
        'containerid': 1076031665372775,
        'page': page
    }
    url = base_url + urlencode(params)
    try:
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            return response.json()
    except requests.ConnectionError as e:
        print('抓取错误', e.args)


# 解析页面返回的json数据
def parse_page(json):
    items = json.get('data').get('cards')
    for item in items:
        item = item.get('mblog')
        if item:
            data = {
                'id': item.get('id'),
                'text': pq(item.get("text")).text(),  # 仅提取内容中的文本
                'attitudes': item.get('attitudes_count'),
                'comments': item.get('comments_count'),
                'reposts': item.get('reposts_count')
            }
            yield data


if __name__ == '__main__':
    for page in range(1, 10):  # 抓取前十页的数据
        json = get_single_page(page)
        results = parse_page(json)
        for result in results:
            print(result)